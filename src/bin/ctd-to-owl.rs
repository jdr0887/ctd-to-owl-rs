#[macro_use]
extern crate log;

use ctd_to_owl_rs::model::*;
use horned_owl::io::owx;
use horned_owl::model::*;
use horned_owl::ontology;
use horned_owl::vocab::WithIRI;
use humantime::format_duration;
use itertools::Itertools;
use rayon::prelude::*;
use std::collections;
use std::error;
use std::fs;
use std::io;
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

    #[structopt(short = "c", long = "chebi-to-mesh", long_help = "chebi to mesh tsv file", required = true, parse(from_os_str))]
    chebi_to_mesh: path::PathBuf,
}
fn main() -> Result<(), Box<dyn error::Error>> {
    let start = time::Instant::now();
    env_logger::init();
    let options = Options::from_args();
    debug!("{:?}", options);

    let chebi_to_mesh_map: collections::HashMap<String, String> = fs::read_to_string(&options.chebi_to_mesh)?
        .lines()
        .map(|a| a.split("\t").map(str::to_owned).collect_vec())
        .map(|vec| {
            assert_eq!(vec.len(), 2);
            (vec[1].to_string(), vec[0].to_string())
        })
        .collect();

    let model = ctd_input_to_model(&options.input)?;

    let mut prefix_mapping = curie::PrefixMapping::default();
    prefix_mapping.add_prefix("owl", "http://www.w3.org/2002/07/owl#").unwrap();
    prefix_mapping.add_prefix("rdf", "http://www.w3.org/1999/02/22-rdf-syntax-ns#").unwrap();
    prefix_mapping.add_prefix("xml", "http://www.w3.org/XML/1998/namespace").unwrap();
    prefix_mapping.add_prefix("xsd", "http://www.w3.org/2001/XMLSchema#").unwrap();
    prefix_mapping.add_prefix("rdfs", "http://www.w3.org/2000/01/rdf-schema#").unwrap();
    prefix_mapping.add_prefix("CHEBI", ctd_to_owl_rs::CHEBI).unwrap();
    prefix_mapping.add_prefix("PMID", ctd_to_owl_rs::PMID).unwrap();
    prefix_mapping.add_prefix("MESH", ctd_to_owl_rs::MESH).unwrap();
    prefix_mapping.add_prefix("NCBITaxon", ctd_to_owl_rs::NCBI_TAXON).unwrap();
    prefix_mapping.add_prefix("NCBIGENE", ctd_to_owl_rs::NCBIGENE).unwrap();
    prefix_mapping.add_prefix("DC", ctd_to_owl_rs::DC).unwrap();

    let output_dir: path::PathBuf = options.output;
    fs::create_dir_all(&output_dir)?;

    model.par_chunks(40000).enumerate().for_each(|(idx, model_chunk)| {
        let ontology = build_ontology(model_chunk.to_vec(), &chebi_to_mesh_map).unwrap();
        let output_path = output_dir.clone().join(format!("{}.owx", idx));
        let output = fs::File::create(&output_path).unwrap();
        info!("writing: {:?}", output_path);
        let mut buf_writer = io::BufWriter::new(output);
        owx::writer::write(&mut buf_writer, &ontology, Some(&prefix_mapping)).unwrap();
    });

    // let ontology = build_ontology(model, &chebi_to_mesh_map).unwrap();
    // let output = fs::File::create(&options.output).unwrap();
    // info!("writing: {:?}", &options.output);
    // let mut buf_writer = io::BufWriter::new(output);
    // owx::writer::write(&mut buf_writer, &ontology, Some(&prefix_mapping)).unwrap();

    info!("Duration: {}", format_duration(start.elapsed()).to_string());
    Ok(())
}

fn build_ontology(model: Vec<IXN>, chebi_to_mesh_map: &collections::HashMap<String, String>) -> Result<ontology::axiom_mapped::AxiomMappedOntology, Box<dyn error::Error>> {
    let build = horned_owl::model::Build::new();
    let mut ontology = ontology::axiom_mapped::AxiomMappedOntology::default();
    let provided_by_ap = build.annotation_property("http://purl.org/pav/providedBy");
    let ontology_root_iri = build.iri("http://ctdbase.org");
    ontology.insert(Axiom::DeclareAnnotationProperty(DeclareAnnotationProperty { 0: provided_by_ap.clone() }));
    ontology.insert(Axiom::OntologyAnnotation(OntologyAnnotation { 0: Annotation { ap: provided_by_ap.clone(), av: AnnotationValue::IRI(ontology_root_iri.clone()) } }));

    for ixn in model.iter() {
        let graph_iri = build.iri(format!("{}{}", ctd_to_owl_rs::CTDIXN, ixn.id));
        ontology.insert(Axiom::AnnotationAssertion(AnnotationAssertion::new(
            graph_iri.clone(),
            Annotation { ap: provided_by_ap.clone(), av: AnnotationValue::IRI(ontology_root_iri.clone()) },
        )));

        for (taxon_idx, taxon) in ixn.taxon.iter().enumerate() {
            let ixn_individual_iri = build.iri(format!("{}#{}", graph_iri.to_string(), taxon_idx));
            debug!("ixn_individual_iri: {:?}", ixn_individual_iri);
            ontology.insert(Axiom::DeclareNamedIndividual(DeclareNamedIndividual { 0: ixn_individual_iri.clone().into() }));

            match process_actor(&build, ixn, &taxon_idx, &taxon, &ixn_individual_iri, chebi_to_mesh_map, &ixn.axns, &ixn.actors) {
                Some((_, actor_axioms)) => {
                    debug!("using ixn: {}", ixn.id);
                    actor_axioms.into_iter().for_each(|axiom| {
                        ontology.insert(axiom);
                    });
                }
                _ => {
                    debug!("skipping ixn: {}", ixn.id);
                }
            }
        }
    }
    Ok(ontology)
}

fn process_actor(
    build: &Build,
    ixn: &IXN,
    taxon_idx: &usize,
    taxon: &Taxon,
    ixn_individual_iri: &IRI,
    chebi_to_mesh_map: &collections::HashMap<String, String>,
    axns: &Vec<AXN>,
    actors: &Vec<Actor>,
) -> Option<(NamedIndividual, Vec<Axiom>)> {
    let codes = axns.iter().map(|a| a.code.clone()).collect_vec();

    if codes.iter().all(|p| p.as_str() == "w" && actors.iter().all(|a| a.actor_type.as_str() != "ixn")) {
        // cotreatment
        let mut axioms: Vec<Axiom> = Vec::new();
        actors.iter().for_each(|actor| {
            let (actor_individual, mut atomic_actor_axioms) =
                get_local_individual_and_axioms(build, actor, taxon_idx, chebi_to_mesh_map).expect("could not get actor class and entity");
            axioms.append(&mut atomic_actor_axioms);
            axioms.append(&mut build_object_property_assertion(&build.object_property(ctd_to_owl_rs::HAS_INPUT.clone()), ixn_individual_iri, &actor_individual).unwrap());
            axioms.append(&mut build_class_assertion(&build.class(ctd_to_owl_rs::COTREATMENT.clone()), ixn_individual_iri).unwrap());
            let mut remnant_axioms = add_remnants(build, ixn, taxon, &ixn_individual_iri).unwrap();
            axioms.append(&mut remnant_axioms);
        });
        return Some((ixn_individual_iri.clone().into(), axioms));
    }

    if codes.iter().all(|p| p.as_str() == "b" && actors.iter().all(|a| a.actor_type.as_str() != "ixn")) {
        // binding
        let mut axioms: Vec<Axiom> = Vec::new();
        actors.iter().for_each(|actor| {
            let (actor_individual, mut atomic_actor_axioms) =
                get_local_individual_and_axioms(build, actor, taxon_idx, chebi_to_mesh_map).expect("could not get actor class and entity");
            axioms.append(&mut atomic_actor_axioms);
            axioms.append(&mut build_object_property_assertion(&build.object_property(ctd_to_owl_rs::HAS_INPUT.clone()), ixn_individual_iri, &actor_individual).unwrap());
            axioms.append(&mut build_class_assertion(&build.class(ctd_to_owl_rs::BINDING.clone()), ixn_individual_iri).unwrap());
            let mut remnant_axioms = add_remnants(build, ixn, taxon, &ixn_individual_iri).unwrap();
            axioms.append(&mut remnant_axioms);
        });
        return Some((ixn_individual_iri.clone().into(), axioms));
    }

    if codes.iter().all(|p| p.as_str() == "rxn") && actors.len() == 2 && actors[1].actor_type == "ixn" {
        let mut axioms: Vec<Axiom> = Vec::new();

        let subject = &actors[0];
        let results = match subject.actor_type.as_str() {
            "ixn" => match process_actor(build, ixn, taxon_idx, taxon, ixn_individual_iri, chebi_to_mesh_map, &subject.axns, &subject.actors) {
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
                match process_actor(build, ixn, taxon_idx, taxon, ixn_individual_iri, chebi_to_mesh_map, &target.axns, &target.actors) {
                    Some((_, target_axioms)) => {
                        debug!("rxn - target.id: {:?}", target.id);

                        axioms.append(&mut target_axioms.clone());
                        axioms.append(&mut subject_axioms.clone());
                        let subject_process_iri = build.iri(format!("{}-process", subject_individual.0));
                        axioms.push(Axiom::DeclareNamedIndividual(DeclareNamedIndividual { 0: subject_process_iri.clone().into() }));

                        let (target_individual, target_axioms) = process_actor(build, ixn, taxon_idx, taxon, ixn_individual_iri, chebi_to_mesh_map, &target.axns, &target.actors)
                            .expect(format!("failed to process actor: {:?}", target).as_str());
                        axioms.append(&mut target_axioms.clone());

                        let process_class = build.class(ctd_to_owl_rs::PROCESS.clone());
                        axioms.push(Axiom::DeclareClass(DeclareClass { 0: process_class.clone() }));
                        axioms.push(Axiom::ClassAssertion(ClassAssertion { ce: process_class.clone().into(), i: ixn_individual_iri.clone().into() }));

                        let part_of_class = build.object_property(ctd_to_owl_rs::PART_OF.clone());
                        axioms.push(Axiom::DeclareObjectProperty(DeclareObjectProperty { 0: part_of_class.clone() }));
                        axioms.push(Axiom::ObjectPropertyAssertion(ObjectPropertyAssertion::new(
                            part_of_class.into(),
                            target_individual.clone(),
                            ixn_individual_iri.clone().into(),
                        )));

                        let axn = ixn.axns.iter().next().expect("could not get AXN from IXN");
                        let process_to_process_op = process_to_process(&build, &axn.degree_code);
                        axioms.push(Axiom::DeclareObjectProperty(DeclareObjectProperty { 0: process_to_process_op.clone() }));
                        axioms.push(Axiom::ObjectPropertyAssertion(ObjectPropertyAssertion::new(
                            process_to_process_op.into(),
                            subject_process_iri.clone().into(),
                            target_individual.clone(),
                        )));
                    }
                    None => {
                        debug!("failed to process actor: {:?}", target)
                    }
                };
            }
            None => {}
        }

        return Some((ixn_individual_iri.clone().into(), axioms));
    }

    if codes.iter().any(|p| ctd_to_owl_rs::AXN_CODES.contains(&p.as_str())) && actors.len() == 2 && actors[1].actor_type != "ixn" {
        let mut axioms: Vec<Axiom> = Vec::new();

        let subject = &actors[0];
        let results = match subject.actor_type.as_str() {
            "ixn" => match process_actor(build, ixn, taxon_idx, taxon, ixn_individual_iri, chebi_to_mesh_map, &subject.axns, &subject.actors) {
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
                axioms.push(Axiom::DeclareNamedIndividual(DeclareNamedIndividual { 0: subject_individual.clone() }));
                axioms.push(Axiom::DeclareNamedIndividual(DeclareNamedIndividual { 0: subject_process.clone() }));

                axioms.append(&mut build_class_assertion(&build.class(ctd_to_owl_rs::PROCESS.clone()), &subject_process.0).unwrap());
                axioms
                    .append(&mut build_object_property_assertion(&build.object_property(ctd_to_owl_rs::HAS_PARTICIPANT.clone()), &subject_process.0, &subject_individual).unwrap());

                axioms.append(&mut subject_axioms.clone());
                let axn = ixn.axns.iter().next().expect("could not get AXN from IXN");

                let target = &actors[1];
                let (target_individual, target_axioms) =
                    get_local_individual_and_axioms(build, &target, taxon_idx, chebi_to_mesh_map).expect("could not get actor class and entity");
                axioms.append(&mut target_axioms.clone());
                let class_map = ctd_to_owl_rs::get_class_map();
                codes.iter().filter(|code| class_map.contains_key(code.as_str())).enumerate().for_each(|(idx, code)| {
                    let ixn_type = class_map.get(code).expect(format!("class not found for code: {:?}", code).as_str());
                    let ixn_type_class = build.class(ixn_type);
                    let local_ixn_iri = build.iri(format!("{}{}#{}-target-{}", ctd_to_owl_rs::CTDIXN, ixn.id, taxon_idx, idx));

                    axioms.append(&mut build_class_assertion(&ixn_type_class, &local_ixn_iri).unwrap());

                    axioms.push(Axiom::DeclareNamedIndividual(DeclareNamedIndividual { 0: local_ixn_iri.clone().into() }));

                    axioms.push(Axiom::ObjectPropertyAssertion(ObjectPropertyAssertion::new(
                        build.object_property(ctd_to_owl_rs::HAS_PARTICIPANT.clone()).into(),
                        local_ixn_iri.clone().into(),
                        target_individual.clone(),
                    )));

                    let process_to_process_op = process_to_process(&build, &axn.degree_code);
                    axioms.append(&mut build_object_property_assertion(&process_to_process_op, &subject_process.0, &local_ixn_iri.clone().into()).unwrap());

                    let part_of_op = build.object_property(ctd_to_owl_rs::PART_OF.clone());
                    axioms.push(Axiom::DeclareObjectProperty(DeclareObjectProperty { 0: part_of_op.clone() }));
                    axioms.push(Axiom::ObjectPropertyAssertion(ObjectPropertyAssertion::new(part_of_op.into(), local_ixn_iri.clone().into(), ixn_individual_iri.clone().into())));
                });
                let mut remnant_axioms = add_remnants(build, ixn, taxon, &ixn_individual_iri).unwrap();
                axioms.append(&mut remnant_axioms);
            }
            None => {}
        }
        return Some((ixn_individual_iri.clone().into(), axioms));
    }

    debug!("not using ixn: {:?}", ixn.id);
    None
}

fn get_local_individual_and_axioms(
    build: &Build,
    actor: &Actor,
    taxon_idx: &usize,
    chebi_to_mesh_map: &collections::HashMap<String, String>,
) -> Result<(NamedIndividual, Vec<Axiom>), Box<dyn error::Error>> {
    let (actor_class, actor_entity, actor_text, actor_label) = match actor.actor_type.as_str() {
        "chemical" => {
            let chebi_mapping = chebi_to_mesh_map.get(actor.id.as_str());
            let actor_class = match chebi_mapping {
                Some(c) => build.class(c.replace("CHEBI:", ctd_to_owl_rs::CHEBI)),
                None => {
                    debug!("no mapping for: {:?}", actor.id);
                    build.class(actor.id.replace("MESH:", format!("{}", ctd_to_owl_rs::MESH).as_str()))
                }
            };

            let actor_text = match &actor.text {
                Some(t) => t.clone(),
                None => {
                    warn!("text is empty - actor: {:?}", actor);
                    String::from("")
                }
            };
            let label = format!("{}#{}-{}", actor_text, actor.parent_id, actor.position);
            (actor_class, build.class(ctd_to_owl_rs::CHEMICAL_ENTITY.clone()), actor_text.clone(), label)
        }
        "gene" => {
            let actor_class = build.class(actor.id.replace("GENE:", ctd_to_owl_rs::NCBIGENE));
            let (actor_text, label) = match &actor.text {
                Some(t) => (t.clone(), format!("{}#{}-{}", t, actor.parent_id, actor.position)),
                None => {
                    let seq_id_value = match &actor.seq_id {
                        Some(s) => s.clone(),
                        None => {
                            warn!("seq_id is empty - actor: {:?}", actor);
                            String::from("")
                        }
                    };
                    (seq_id_value.clone(), format!("{}#{}-{}", seq_id_value, actor.parent_id, actor.position))
                }
            };
            (actor_class, build.class(ctd_to_owl_rs::GENE_ENTITY.clone()), actor_text, label)
        }
        _ => {
            panic!("should never get here")
        }
    };

    let mut axioms: Vec<Axiom> = Vec::new();

    let actor_individual_iri = build.iri(format!("{}{}#{}-{}", ctd_to_owl_rs::CTDIXN, actor.parent_id, taxon_idx, actor.position));
    axioms.push(Axiom::DeclareNamedIndividual(DeclareNamedIndividual { 0: actor_individual_iri.clone().into() }));

    // actorInd Type actorClass,
    axioms.append(&mut build_class_assertion(&actor_class, &actor_individual_iri)?);

    // actorClass Annotation(RDFSLabel, typeLabel),
    axioms.push(Axiom::AnnotationAssertion(AnnotationAssertion::new(
        actor_class.0.clone(),
        Annotation { ap: build.annotation_property(horned_owl::vocab::RDFS::Label.iri_s()), av: AnnotationValue::Literal(Literal::Simple { literal: actor_text.clone() }) },
    )));

    // actorInd Type nodeType,
    axioms.append(&mut build_class_assertion(&actor_entity, &actor_individual_iri)?);

    // actorInd Annotation(RDFSLabel, label)
    axioms.push(Axiom::AnnotationAssertion(AnnotationAssertion::new(
        actor_individual_iri.clone(),
        Annotation { ap: build.annotation_property(horned_owl::vocab::RDFS::Label.iri_s()), av: AnnotationValue::Literal(Literal::Simple { literal: actor_label }) },
    )));

    match &actor.form {
        Some(s) => {
            let form_ap = build.annotation_property("http://ctd.example.org/has_form");
            axioms.push(Axiom::DeclareAnnotationProperty(DeclareAnnotationProperty { 0: form_ap.clone() }));
            axioms.push(Axiom::AnnotationAssertion(AnnotationAssertion::new(
                actor_individual_iri.clone(),
                Annotation { ap: form_ap.clone(), av: AnnotationValue::Literal(Literal::Simple { literal: s.to_string() }) },
            )));
        }
        _ => {}
    };
    Ok((actor_individual_iri.into(), axioms))
}

fn add_remnants(build: &Build, ixn: &IXN, taxon: &Taxon, ixn_individual_iri: &IRI) -> Result<Vec<Axiom>, Box<dyn error::Error>> {
    let mut axioms: Vec<Axiom> = Vec::new();
    let pm_ids = ixn.reference.iter().map(|r| format!("{}/{}", ctd_to_owl_rs::PMID, r.pm_id)).collect_vec();

    axioms.push(Axiom::DeclareAnnotationProperty(DeclareAnnotationProperty { 0: build.annotation_property(ctd_to_owl_rs::DC_SOURCE.clone()) }));

    pm_ids.iter().for_each(|pm_id_iri| {
        axioms.push(Axiom::AnnotationAssertion(AnnotationAssertion::new(
            ixn_individual_iri.clone(),
            Annotation { ap: build.annotation_property(ctd_to_owl_rs::DC_SOURCE.clone()), av: AnnotationValue::IRI(build.iri(pm_id_iri.clone())) },
        )))
    });
    let organism_iri = build.iri(format!("{}-organism", ixn_individual_iri.to_string()));
    let taxon_iri = build.iri(format!("{}{}", ctd_to_owl_rs::NCBI_TAXON, &taxon.id));
    axioms.append(&mut build_class_assertion(&taxon_iri.into(), &organism_iri)?);
    axioms.push(Axiom::AnnotationAssertion(AnnotationAssertion::new(
        organism_iri.clone(),
        Annotation {
            ap: build.annotation_property(horned_owl::vocab::RDFS::Label.iri_s()),
            av: AnnotationValue::Literal(Literal::Simple { literal: format!("{}#{}", &taxon.text, &ixn.id) }),
        },
    )));

    axioms.push(Axiom::DeclareObjectProperty(DeclareObjectProperty { 0: build.object_property(ctd_to_owl_rs::OCCURS_IN.clone()) }));
    axioms.push(Axiom::DeclareNamedIndividual(DeclareNamedIndividual { 0: organism_iri.clone().into() }));
    axioms.push(Axiom::ObjectPropertyAssertion(ObjectPropertyAssertion::new(
        ObjectPropertyExpression::ObjectProperty(build.object_property(ctd_to_owl_rs::OCCURS_IN.clone())),
        ixn_individual_iri.clone().into(),
        organism_iri.clone().into(),
    )));
    Ok(axioms)
}

fn ctd_input_to_model(ctd_input_path: &path::PathBuf) -> Result<Vec<IXN>, Box<dyn error::Error>> {
    let data = fs::read_to_string(ctd_input_path)?;
    let ixnset_element = Element::parse(data.as_bytes())?;
    parse_input(&ixnset_element)
}

fn build_object_property_assertion(object_property: &ObjectProperty, ixn_individual_iri: &IRI, actor_individual: &NamedIndividual) -> Result<Vec<Axiom>, Box<dyn error::Error>> {
    let mut axioms = Vec::new();
    axioms.push(Axiom::DeclareObjectProperty(DeclareObjectProperty { 0: object_property.clone() }));
    axioms.push(Axiom::ObjectPropertyAssertion(ObjectPropertyAssertion::new(object_property.clone().into(), ixn_individual_iri.clone().into(), actor_individual.clone())));
    Ok(axioms)
}

fn build_class_assertion(class: &Class, ixn_individual_iri: &IRI) -> Result<Vec<Axiom>, Box<dyn error::Error>> {
    let mut axioms = Vec::new();
    axioms.push(Axiom::DeclareClass(DeclareClass { 0: class.clone() }));
    axioms.push(Axiom::ClassAssertion(ClassAssertion { ce: class.clone().into(), i: ixn_individual_iri.clone().into() }));
    Ok(axioms)
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
                                parse_actor_element(&mut ixn.actors, &ixn_child_element);
                            }
                            _ => {}
                        },
                        _ => {}
                    }
                }
                debug!("{:?}", ixn);
                ixns.push(ixn);
            }
            _ => {}
        }
    }
    Ok(ixns)
}

fn parse_actor_element(actors: &mut Vec<Actor>, element: &Element) {
    let mut actor = get_actor_from_element(element).unwrap();
    match element.get_text() {
        Some(s) => {
            actor.text = Some(s.to_string());
            actors.push(actor);
        }
        None => {
            if !element.children.is_empty() {
                for a_node in element.children.iter() {
                    match a_node {
                        XMLNode::Element(a_element) => match a_element.name.as_str() {
                            "axn" => {
                                let a_created_axn = get_axn_from_element(a_element).unwrap();
                                actor.axns.push(a_created_axn.clone());
                            }
                            "actor" => parse_actor_element(&mut actor.actors, a_element),
                            _ => {}
                        },
                        _ => {}
                    }
                }
            }
            actors.push(actor);
        }
    }
}
