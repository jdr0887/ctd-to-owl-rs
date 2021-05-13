#[macro_use]
extern crate log;

use ctd_to_owl_rs::model::*;
use ctd_to_owl_rs::OBO;
use horned_owl::io::owx;
use horned_owl::model::*;
use horned_owl::ontology;
use horned_owl::vocab::WithIRI;
use humantime::format_duration;
use itertools::{all, Itertools};
use std::collections;
use std::collections::{HashMap, HashSet};
use std::error;
use std::fs;
use std::io;
use std::io::BufRead;
use std::path;
use std::time;
use structopt::StructOpt;
use xmltree::{Element, XMLNode};

#[derive(StructOpt, Debug)]
#[structopt(name = "ctd-to-owl-rs", about = "convert ctd xml to owx")]
struct Options {
    #[structopt(short = "i", long = "input", long_help = "input", required = true, parse(from_os_str))]
    input: path::PathBuf,

    #[structopt(short = "o", long = "output", long_help = "output", required = true, parse(from_os_str))]
    output: path::PathBuf,
}
fn main() -> Result<(), Box<dyn error::Error>> {
    let start = time::Instant::now();
    env_logger::init();
    let options = Options::from_args();
    debug!("{:?}", options);

    let chebi_to_mesh_map: collections::HashMap<String, String> = include_str!("../data/chebi_mesh.tsv")
        .lines()
        .map(|a| a.split("\t").map(str::to_owned).collect_vec())
        .map(|vec| {
            assert_eq!(vec.len(), 2);
            (vec[1].to_string(), vec[0].to_string())
        })
        .collect();

    let build = horned_owl::model::Build::new();

    let data = fs::read_to_string(&options.input)?;
    let ixnset_element = Element::parse(data.as_bytes()).unwrap();

    let model = parse_input(&ixnset_element)?;

    let mut ontology = ontology::axiom_mapped::AxiomMappedOntology::default();

    ontology.insert(Axiom::OntologyAnnotation(OntologyAnnotation {
        0: Annotation { ap: build.annotation_property("http://purl.org/pav/providedBy"), av: AnnotationValue::IRI(build.iri("http://ctdbase.org")) },
    }));

    for ixn in model.iter() {
        // let graph_iri = format!("{}{}", "http://ctdbase.org/detail.go?type=relationship&ixnId=", ixn.id);
        for (taxon_idx, taxon) in ixn.taxon.iter().enumerate() {
            match process_actor(&build, ixn, &taxon_idx, taxon, &chebi_to_mesh_map, &ixn.axns, &ixn.actors) {
                Some((actor_individual, actor_axioms)) => {
                    info!("using ixn: {}", ixn.id);
                    actor_axioms.into_iter().for_each(|axiom| {
                        ontology.insert(axiom);
                    });
                }
                _ => {
                    warn!("skipping ixn: {}", ixn.id);
                }
            }
        }
    }

    let mut prefix_mapping = curie::PrefixMapping::default();
    prefix_mapping.add_prefix("owl", "http://www.w3.org/2002/07/owl#").unwrap();
    prefix_mapping.add_prefix("rdf", "http://www.w3.org/1999/02/22-rdf-syntax-ns#").unwrap();
    prefix_mapping.add_prefix("xml", "http://www.w3.org/XML/1998/namespace").unwrap();
    prefix_mapping.add_prefix("xsd", "http://www.w3.org/2001/XMLSchema#").unwrap();
    prefix_mapping.add_prefix("rdfs", "http://www.w3.org/2000/01/rdf-schema#").unwrap();
    prefix_mapping.add_prefix("CHEBI", ctd_to_owl_rs::CHEBI).unwrap();
    // prefix_mapping.add_prefix("PMID", ctd_to_owl_rs::PMID).unwrap();
    // prefix_mapping.add_prefix("MESH", ctd_to_owl_rs::MESH).unwrap();
    prefix_mapping.add_prefix("NCBITaxon", ctd_to_owl_rs::NCBI_TAXON).unwrap();
    prefix_mapping.add_prefix("NCBIGENE", ctd_to_owl_rs::NCBIGENE).unwrap();

    let output = fs::File::create(&options.output).ok().unwrap();
    info!("writing: {:?}", &options.output);
    let mut buf_writer = io::BufWriter::new(output);
    owx::writer::write(&mut buf_writer, &ontology, Some(&prefix_mapping))?;

    info!("Duration: {}", format_duration(start.elapsed()).to_string());
    Ok(())
}

fn process_actor(
    build: &Build,
    ixn: &IXN,
    taxon_idx: &usize,
    taxon: &Taxon,
    chebi_to_mesh_map: &collections::HashMap<String, String>,
    axns: &Vec<AXN>,
    actors: &Vec<Actor>,
) -> Option<(NamedIndividual, Vec<Axiom>)> {
    debug!("actors: {:?}", actors);
    debug!("axns: {:?}", axns);

    let codes = axns.iter().map(|a| a.code.clone()).collect_vec();
    let ixn_individual = build.iri(format!("{}{}#{}", ctd_to_owl_rs::CTDIXN, ixn.id, taxon_idx));

    if codes.iter().all(|p| p.as_str() == "w" && actors.iter().all(|a| a.actor_type.as_str() != "ixn")) {
        // cotreatment
        let mut axioms: Vec<Axiom> = Vec::new();
        actors.iter().for_each(|actor| {
            let (actor_individual, mut atomic_actor_axioms) =
                get_local_individual_and_axioms(build, actor, taxon_idx, chebi_to_mesh_map).expect("could not get actor class and entity");
            axioms.append(&mut atomic_actor_axioms);
            axioms.push(Axiom::ObjectPropertyAssertion(ObjectPropertyAssertion::new(
                build.object_property(ctd_to_owl_rs::HAS_INPUT.clone()).into(),
                NamedIndividual::from(ixn_individual.clone()),
                actor_individual.clone(),
            )));
            axioms.push(Axiom::ClassAssertion(ClassAssertion {
                ce: ClassExpression::from(build.class(ctd_to_owl_rs::COTREATMENT.clone())),
                i: NamedIndividual::from(ixn_individual.clone()),
            }));
            let mut remnant_axioms = add_remnants(build, ixn, &ixn_individual, taxon).unwrap();
            axioms.append(&mut remnant_axioms);
        });
        return Some((NamedIndividual::from(ixn_individual.clone()), axioms));
    }

    if codes.iter().all(|p| p.as_str() == "b" && actors.iter().all(|a| a.actor_type.as_str() != "ixn")) {
        // binding
        let mut axioms: Vec<Axiom> = Vec::new();
        actors.iter().for_each(|actor| {
            let (actor_individual, mut atomic_actor_axioms) =
                get_local_individual_and_axioms(build, actor, taxon_idx, chebi_to_mesh_map).expect("could not get actor class and entity");
            axioms.append(&mut atomic_actor_axioms);
            axioms.push(Axiom::ObjectPropertyAssertion(ObjectPropertyAssertion::new(
                build.object_property(ctd_to_owl_rs::HAS_INPUT.clone()).into(),
                NamedIndividual::from(ixn_individual.clone()),
                actor_individual.clone(),
            )));
            axioms.push(Axiom::ClassAssertion(ClassAssertion {
                ce: ClassExpression::from(build.class(ctd_to_owl_rs::BINDING.clone())),
                i: NamedIndividual::from(ixn_individual.clone()),
            }));
            let mut remnant_axioms = add_remnants(build, ixn, &ixn_individual, taxon).unwrap();
            axioms.append(&mut remnant_axioms);
        });
        return Some((NamedIndividual::from(ixn_individual.clone()), axioms));
    }

    if codes.iter().all(|p| p.as_str() == "rxn") && actors.len() == 2 && actors[1].actor_type == "ixn" {
        let mut axioms: Vec<Axiom> = Vec::new();
        let subject = &actors[0];
        let results = match subject.actor_type.as_str() {
            "ixn" => match process_actor(build, ixn, taxon_idx, taxon, chebi_to_mesh_map, &subject.axns, &subject.actors) {
                Some((subject_individual, subject_axioms)) => Some((subject_individual, subject_axioms)),
                None => None,
            },
            _ => {
                let (subject_individual, subject_axioms) =
                    get_local_individual_and_axioms(build, subject, taxon_idx, chebi_to_mesh_map).expect("could not get actor class and entity");
                Some((subject_individual, subject_axioms))
            }
        };

        match results {
            Some((subject_individual, subject_axioms)) => {
                let target = &actors[1];
                match process_actor(build, ixn, taxon_idx, taxon, chebi_to_mesh_map, &target.axns, &target.actors) {
                    Some((target_individual, target_axioms)) => {
                        debug!("rxn - target.id: {:?}", target.id);

                        axioms.append(&mut target_axioms.clone());
                        axioms.append(&mut subject_axioms.clone());
                        let subject_process = build.named_individual(format!("{}-process", subject_individual.0));
                        let (target_individual, target_axioms) = process_actor(build, ixn, taxon_idx, taxon, chebi_to_mesh_map, &target.axns, &target.actors)
                            .expect(format!("failed to process actor: {:?}", target).as_str());
                        axioms.append(&mut target_axioms.clone());
                        axioms.push(Axiom::ClassAssertion(ClassAssertion {
                            ce: ClassExpression::from(build.class(ctd_to_owl_rs::PROCESS.clone())),
                            i: NamedIndividual::from(ixn_individual.clone()),
                        }));

                        axioms.push(Axiom::ObjectPropertyAssertion(ObjectPropertyAssertion::new(
                            build.object_property(ctd_to_owl_rs::PART_OF.clone()).into(),
                            target_individual.clone(),
                            NamedIndividual::from(ixn_individual.clone()),
                        )));

                        let axn = ixn.axns.iter().next().expect("could not get AXN from IXN");
                        axioms.push(Axiom::ObjectPropertyAssertion(ObjectPropertyAssertion::new(
                            process_to_process(&build, &axn.degree_code).into(),
                            subject_process.clone(),
                            target_individual.clone(),
                        )));
                    }
                    None => {
                        warn!("failed to process actor: {:?}", target)
                    }
                };
            }
            None => {}
        }

        return Some((NamedIndividual::from(ixn_individual.clone()), axioms));
    }

    if codes.iter().any(|p| ctd_to_owl_rs::AXN_CODES.contains(&p.as_str())) && actors.len() == 2 && actors[1].actor_type != "ixn" {
        let mut axioms: Vec<Axiom> = Vec::new();

        let subject = &actors[0];
        let results = match subject.actor_type.as_str() {
            "ixn" => match process_actor(build, ixn, taxon_idx, taxon, chebi_to_mesh_map, &subject.axns, &subject.actors) {
                Some((subject_individual, subject_axioms)) => Some((subject_individual, subject_axioms)),
                None => None,
            },
            _ => {
                let (subject_individual, subject_axioms) =
                    get_local_individual_and_axioms(build, subject, taxon_idx, chebi_to_mesh_map).expect("could not get actor class and entity");
                Some((subject_individual, subject_axioms))
            }
        };

        match results {
            Some((subject_individual, subject_axioms)) => {
                let subject_process = build.named_individual(format!("{}-process", subject_individual.0));

                axioms.push(Axiom::ClassAssertion(ClassAssertion { ce: ClassExpression::from(build.class(ctd_to_owl_rs::PROCESS.clone())), i: subject_process.clone() }));
                axioms.push(Axiom::ObjectPropertyAssertion(ObjectPropertyAssertion::new(
                    build.object_property(ctd_to_owl_rs::HAS_PARTICIPANT.clone()).into(),
                    subject_process.clone(),
                    subject_individual.clone(),
                )));
                axioms.append(&mut subject_axioms.clone());
                let axn = ixn.axns.iter().next().expect("could not get AXN from IXN");

                let target = &actors[1];
                let (target_individual, target_axioms) =
                    get_local_individual_and_axioms(build, &target, taxon_idx, chebi_to_mesh_map).expect("could not get actor class and entity");
                axioms.append(&mut target_axioms.clone());
                let class_map = ctd_to_owl_rs::get_class_map();
                codes.iter().filter(|code| class_map.contains_key(code.as_str())).enumerate().for_each(|(idx, code)| {
                    let ixn_type = class_map.get(code).expect(format!("class not found for code: {:?}", code).as_str());
                    let local_ixn_individual = build.named_individual(format!("{}{}#{}-target-{}", ctd_to_owl_rs::CTDIXN, ixn.id, taxon_idx, idx));

                    axioms.push(Axiom::ClassAssertion(ClassAssertion { ce: ClassExpression::from(build.class(ixn_type)), i: local_ixn_individual.clone() }));
                    axioms.push(Axiom::ObjectPropertyAssertion(ObjectPropertyAssertion::new(
                        build.object_property(ctd_to_owl_rs::HAS_PARTICIPANT.clone()).into(),
                        local_ixn_individual.clone(),
                        target_individual.clone(),
                    )));
                    axioms.push(Axiom::ObjectPropertyAssertion(ObjectPropertyAssertion::new(
                        process_to_process(&build, &axn.degree_code).into(),
                        subject_process.clone(),
                        local_ixn_individual.clone(),
                    )));
                    axioms.push(Axiom::ObjectPropertyAssertion(ObjectPropertyAssertion::new(
                        build.object_property(ctd_to_owl_rs::PART_OF.clone()).into(),
                        local_ixn_individual.clone(),
                        NamedIndividual::from(ixn_individual.clone()),
                    )));
                });
                let mut remnant_axioms = add_remnants(build, ixn, &ixn_individual, taxon).unwrap();
                axioms.append(&mut remnant_axioms);
            }
            None => {}
        }
        return Some((NamedIndividual::from(ixn_individual.clone()), axioms));
    }

    warn!("not handling ixn: {:?}", ixn.id);
    None
}

fn get_local_individual_and_axioms(
    build: &Build,
    actor: &Actor,
    taxon_idx: &usize,
    chebi_to_mesh_map: &HashMap<String, String>,
) -> Result<(NamedIndividual, Vec<Axiom>), Box<dyn error::Error>> {
    let (actor_class, actor_entity, actor_text, actor_label) = match actor.actor_type.as_str() {
        "chemical" => {
            let chebi_mapping = chebi_to_mesh_map.get(actor.id.as_str());
            let actor_class = match chebi_mapping {
                Some(c) => build.class(c.replace("CHEBI:", ctd_to_owl_rs::CHEBI)),
                None => {
                    warn!("no mapping for: {:?}", actor.id);
                    build.class(actor.id.replace("MESH:", format!("{}", ctd_to_owl_rs::MESH).as_str()))
                }
            };
            let actor_text = actor.text.as_ref().expect("failed to get actor text");
            let label = format!("{}#{}-{}", actor_text, actor.parent_id, actor.position);
            (actor_class, build.class(format!("{}/{}", ctd_to_owl_rs::OBO, "CHEBI_24431")), actor_text.clone(), label)
        }
        "gene" => {
            let actor_class = build.class(actor.id.replace("GENE:", ctd_to_owl_rs::NCBIGENE));
            let (actor_text, label) = match &actor.text {
                Some(t) => (t, format!("{}#{}-{}", t, actor.parent_id, actor.position)),
                None => {
                    let t = actor.seq_id.as_ref().unwrap();
                    (t, format!("{}#{}-{}", t, actor.parent_id, actor.position))
                }
            };
            (actor_class, build.class(format!("{}/{}", ctd_to_owl_rs::OBO, "SO_0000704")), actor_text.to_string(), label)
        }
        _ => {
            panic!("should never get here")
        }
    };

    let actor_individual_iri = build.iri(format!("{}{}#{}-{}", ctd_to_owl_rs::CTDIXN, actor.parent_id, taxon_idx, actor.position));
    let mut axioms: Vec<Axiom> = Vec::new();

    // actorInd Type actorClass,
    axioms.push(Axiom::ClassAssertion(ClassAssertion { ce: ClassExpression::from(actor_class.clone()), i: NamedIndividual::from(actor_individual_iri.clone()) }));

    // actorClass Annotation(RDFSLabel, typeLabel),
    axioms.push(Axiom::AnnotationAssertion(AnnotationAssertion::new(
        actor_individual_iri.clone(),
        Annotation { ap: build.annotation_property(horned_owl::vocab::RDFS::Label.iri_s()), av: AnnotationValue::Literal(Literal::Simple { literal: actor_text.clone() }) },
    )));

    // actorInd Type nodeType,
    axioms.push(Axiom::ClassAssertion(ClassAssertion { ce: ClassExpression::from(actor_entity.clone()), i: NamedIndividual::from(actor_individual_iri.clone()) }));

    // actorInd Annotation(RDFSLabel, label)
    axioms.push(Axiom::AnnotationAssertion(AnnotationAssertion::new(
        actor_individual_iri.clone(),
        Annotation { ap: build.annotation_property(horned_owl::vocab::RDFS::Label.iri_s()), av: AnnotationValue::Literal(Literal::Simple { literal: actor_label }) },
    )));

    match &actor.form {
        Some(s) => {
            axioms.push(Axiom::AnnotationAssertion(AnnotationAssertion::new(
                actor_individual_iri.clone(),
                Annotation { ap: build.annotation_property("http://ctd.example.org/has_form"), av: AnnotationValue::Literal(Literal::Simple { literal: s.to_string() }) },
            )));
        }
        _ => {}
    };
    Ok((NamedIndividual::from(actor_individual_iri), axioms))
}

fn add_remnants(build: &Build, ixn: &IXN, ixn_individual: &IRI, taxon: &Taxon) -> Result<Vec<Axiom>, Box<dyn error::Error>> {
    let mut axioms: Vec<Axiom> = Vec::new();
    let pm_ids = ixn.reference.iter().map(|r| format!("{}/{}", ctd_to_owl_rs::PMID, r.pm_id)).collect_vec();

    pm_ids.iter().for_each(|pm_id_iri| {
        axioms.push(Axiom::AnnotationAssertion(AnnotationAssertion::new(
            ixn_individual.clone(),
            Annotation { ap: build.annotation_property(ctd_to_owl_rs::DC_SOURCE.clone()), av: AnnotationValue::IRI(build.iri(pm_id_iri.clone())) },
        )))
    });
    let organism = build.iri(format!("{}#organism", ixn_individual.to_string()));
    let taxon_class = build.class(format!("{}{}", ctd_to_owl_rs::NCBI_TAXON, &taxon.id));
    axioms.push(Axiom::ClassAssertion(ClassAssertion { ce: ClassExpression::from(taxon_class), i: NamedIndividual::from(organism.clone()) }));
    axioms.push(Axiom::AnnotationAssertion(AnnotationAssertion::new(
        organism.clone(),
        Annotation {
            ap: build.annotation_property(horned_owl::vocab::RDFS::Label.iri_s()),
            av: AnnotationValue::Literal(Literal::Simple { literal: format!("{}#{}", &taxon.text, &ixn.id) }),
        },
    )));

    axioms.push(Axiom::ObjectPropertyAssertion(ObjectPropertyAssertion::new(
        ObjectPropertyExpression::ObjectProperty(build.object_property(ctd_to_owl_rs::OCCURS_IN.clone())),
        NamedIndividual::from(ixn_individual.clone()),
        NamedIndividual::from(organism),
    )));
    Ok(axioms)
}

fn parse_input(ixnset_element: &Element) -> Result<Vec<IXN>, Box<dyn error::Error>> {
    let mut ixns = Vec::new();

    for ixn_node in ixnset_element.children.iter() {
        match ixn_node {
            XMLNode::Element(ixn_element) => {
                let ixn_id = ixn_element.attributes.get("id").unwrap();
                let mut ixn = IXN::new();
                ixn.id = ixn_id.parse::<i32>().unwrap();

                for ixn_element_child_node in ixn_element.children.iter() {
                    match ixn_element_child_node {
                        XMLNode::Element(ixn_child_element) => match ixn_child_element.name.as_str() {
                            "taxon" => {
                                let taxon_id = ixn_child_element.attributes.get("id").unwrap();
                                let taxon_text = ixn_child_element.get_text().unwrap().to_string();
                                let taxon = Taxon::new(taxon_id.parse::<i32>().unwrap(), taxon_text);
                                ixn.taxon.push(taxon);
                                // debug!("{:?}", ixn);
                            }
                            "reference" => {
                                let reference_pm_id = ixn_child_element.attributes.get("pmid").unwrap();
                                let reference = Reference::new(reference_pm_id.parse::<i32>().unwrap());
                                ixn.reference.push(reference);
                                // debug!("{:?}", ixn);
                            }
                            "axn" => {
                                let axn = get_axn_from_element(ixn_child_element).unwrap();
                                if ixn.id == axn.parent_id {
                                    ixn.axns.push(axn);
                                }
                                // debug!("{:?}", ixn);
                            }
                            "actor" => {
                                let mut actor = get_actor_from_element(ixn_child_element).unwrap();
                                match ixn_child_element.get_text() {
                                    Some(s) => {
                                        actor.text = Some(s.to_string());
                                        if ixn.id == actor.parent_id {
                                            ixn.actors.push(actor);
                                        }
                                        // debug!("{:?}", ixn);
                                    }
                                    None => {
                                        if ixn.id == actor.parent_id {
                                            ixn.actors.push(actor);
                                        }
                                        // debug!("{:?}", ixn);

                                        for a_node in ixn_child_element.children.iter() {
                                            match a_node {
                                                XMLNode::Element(a_element) => match a_element.name.as_str() {
                                                    "axn" => {
                                                        let a_created_axn = get_axn_from_element(a_element).unwrap();
                                                        let found_actor =
                                                            ixn.actors.iter_mut().find(|b_actor| b_actor.id == a_created_axn.parent_id.to_string()).expect("could not find actor");
                                                        found_actor.axns.push(a_created_axn.clone());
                                                        // debug!("{:?}", ixn);
                                                    }
                                                    "actor" => {
                                                        let mut a_created_actor = get_actor_from_element(a_element).unwrap();
                                                        let a_actor_text = a_element.get_text();
                                                        match a_actor_text {
                                                            Some(t) => {
                                                                a_created_actor.text = Some(t.to_string());
                                                                let found_actor = ixn
                                                                    .actors
                                                                    .iter_mut()
                                                                    .find(|b_actor| b_actor.id == a_created_actor.parent_id.to_string())
                                                                    .expect("could not find actor");
                                                                found_actor.actors.push(a_created_actor.clone());
                                                                // debug!("{:?}", ixn);
                                                            }
                                                            None => {
                                                                let found_actor = ixn
                                                                    .actors
                                                                    .iter_mut()
                                                                    .find(|b_actor| b_actor.id == a_created_actor.parent_id.to_string())
                                                                    .expect("could not find actor");
                                                                found_actor.actors.push(a_created_actor.clone());
                                                                // debug!("{:?}", ixn);

                                                                for b_node in a_element.children.iter() {
                                                                    match b_node {
                                                                        XMLNode::Element(b_element) => match b_element.name.as_str() {
                                                                            "axn" => {
                                                                                let b_created_axn = get_axn_from_element(b_element).unwrap();
                                                                                ixn.actors.iter_mut().for_each(|a_actor| {
                                                                                    let found_actor =
                                                                                        a_actor.actors.iter_mut().find(|b_actor| b_actor.id == b_created_axn.parent_id.to_string());
                                                                                    match found_actor {
                                                                                        Some(a) => a.axns.push(b_created_axn.clone()),
                                                                                        _ => {
                                                                                            // debug!("actor not found for: {:?}", b_created_axn);
                                                                                        }
                                                                                    }
                                                                                });
                                                                            }
                                                                            "actor" => {
                                                                                let mut b_created_actor = get_actor_from_element(b_element).unwrap();
                                                                                let b_actor_text = b_element.get_text();
                                                                                match b_actor_text {
                                                                                    Some(t) => {
                                                                                        b_created_actor.text = Some(t.to_string());

                                                                                        ixn.actors.iter_mut().for_each(|a_axn| {
                                                                                            let found_actor = a_axn
                                                                                                .actors
                                                                                                .iter_mut()
                                                                                                .find(|b_actor| b_actor.id == b_created_actor.parent_id.to_string());
                                                                                            match found_actor {
                                                                                                Some(a) => a.actors.push(b_created_actor.clone()),
                                                                                                _ => {
                                                                                                    // debug!("actor not found for: {:?}", b_created_actor);
                                                                                                }
                                                                                            }
                                                                                        });
                                                                                    }
                                                                                    None => {}
                                                                                }
                                                                                // debug!("{:?}", ixn);
                                                                            }
                                                                            _ => {}
                                                                        },
                                                                        _ => {}
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                    _ => {}
                                                },
                                                _ => {}
                                            }
                                        }
                                    }
                                }
                                // debug!("{:?}", ixn);
                            }
                            _ => {}
                        },
                        _ => {}
                    }
                }
                ixns.push(ixn);
            }
            _ => {}
        }
    }
    Ok(ixns)
}

fn get_actor_from_element(element: &Element) -> Result<Actor, Box<dyn error::Error>> {
    let actor_type = element.attributes.get("type").unwrap();
    let actor_id = element.attributes.get("id").unwrap();
    let actor_position = element.attributes.get("position").unwrap();
    let actor_parent_id = element.attributes.get("parentid").unwrap();
    let form = element.attributes.get("form").cloned();
    //let form_qualifier = element.attributes.get("form_qualifier").cloned();
    let seq_id = element.attributes.get("seqid").cloned();
    // let actor =
    //     Actor::new(actor_type.to_string(), actor_id.to_string(), actor_position.parse::<i32>().unwrap(), actor_parent_id.parse::<i32>().unwrap(), form, form_qualifier, seq_id);
    let actor = Actor::new(actor_type.to_string(), actor_id.to_string(), actor_position.parse::<i8>().unwrap(), actor_parent_id.parse::<i32>().unwrap(), form, None, seq_id);
    Ok(actor)
}

fn get_axn_from_element(element: &Element) -> Result<AXN, Box<dyn error::Error>> {
    let axn_code = element.attributes.get("code").unwrap();
    let axn_degreecode = element.attributes.get("degreecode").unwrap();
    let axn_position = element.attributes.get("position").unwrap();
    let axn_parent_id = element.attributes.get("parentid").unwrap();
    let axn_text = element.get_text().unwrap().to_string();
    let axn = AXN::new(axn_code.into(), axn_degreecode.chars().next().unwrap(), axn_position.parse::<i8>().unwrap(), axn_parent_id.parse::<i32>().unwrap(), axn_text);
    Ok(axn)
}

fn process_to_process(build: &horned_owl::model::Build, degree: &char) -> horned_owl::model::ObjectProperty {
    match degree {
        '1' => build.object_property(ctd_to_owl_rs::CAUSALLY_UPSTREAM_OF.clone()),
        //"0" => // has no effect? never used in CTD_chem_gene_ixns_structured.xml
        '+' => build.object_property(ctd_to_owl_rs::CAUSALLY_UPSTREAM_OF_POSITIVE_EFFECT.clone()),
        '-' => build.object_property(ctd_to_owl_rs::CAUSALLY_UPSTREAM_OF_NEGATIVE_EFFECT.clone()),
        _ => {
            panic!("invalid degree")
        }
    }
}
