#[macro_use]
extern crate log;
//extern crate quick_xml;
extern crate serde_xml_rs;

use ctd_to_owl_rs::model::*;
use ctd_to_owl_rs::{get_class_map, OBO};
use horned_owl::io::owx;
use horned_owl::model::*;
use horned_owl::ontology;
use horned_owl::vocab::WithIRI;
use humantime::format_duration;
use itertools::{all, Itertools};
use std::collections;
use std::error;
use std::fs;
use std::io;
use std::io::BufRead;
use std::path;
use std::time;
use structopt::StructOpt;
use xmltree::{Element, XMLNode};

#[derive(StructOpt, Debug)]
#[structopt(name = "ctd-to-owl-rs", about = "asdfasdf")]
struct Options {
    #[structopt(short = "i", long = "input", long_help = "input", required = true, parse(from_os_str))]
    input: path::PathBuf,

    #[structopt(short = "o", long = "output", long_help = "output", required = true, parse(from_os_str))]
    output: path::PathBuf,

    #[structopt(short = "m", long = "chebi_to_mesh", long_help = "chebi to mesh", required = true, parse(from_os_str))]
    chebi_to_mesh: path::PathBuf,
}
fn main() -> Result<(), Box<dyn error::Error>> {
    let start = time::Instant::now();
    env_logger::init();
    let options = Options::from_args();
    debug!("{:?}", options);

    let chebi_to_mesh_map: collections::HashMap<String, String> = io::BufReader::new(fs::File::open(&options.chebi_to_mesh)?)
        .lines()
        .map(|line| line.unwrap())
        .map(|a| a.split("\t").map(str::to_owned).collect_vec())
        .map(|vec| {
            assert_eq!(vec.len(), 2);
            (vec[0].to_string(), vec[1].to_string())
        })
        .collect();

    let build = horned_owl::model::Build::new();

    let acts_upstream_of = build.object_property(format!("{}/{}", OBO, "RO_0002263"));
    let acts_upstream_of_positive_effect = build.object_property(format!("{}/{}", OBO, "RO_0004034"));
    let acts_upstream_of_negative_effect = build.object_property(format!("{}/{}", OBO, "RO_0004035"));
    let causally_upstream_of = build.object_property(format!("{}/{}", OBO, "RO_0002411"));
    let causally_upstream_of_positive_effect = build.object_property(format!("{}/{}", OBO, "RO_0002304"));
    let causally_upstream_of_negative_effect = build.object_property(format!("{}/{}", OBO, "RO_0002305"));
    let part_of = build.object_property(format!("{}/{}", OBO, "BFO_0000050"));
    let has_participant = build.object_property(format!("{}/{}", OBO, "RO_0000057"));
    let input_of = build.object_property(format!("{}/{}", OBO, "RO_0002352"));
    let enables = build.object_property(format!("{}/{}", OBO, "RO_0002327"));
    let enabled_by = build.object_property(format!("{}/{}", OBO, "RO_0002333"));
    let transports_or_maintains_localization_of = build.object_property(format!("{}/{}", OBO, "RO_0002313"));

    let provided_by = build.annotation_property("http://purl.org/pav/providedBy");

    let class_map = get_class_map(&build)?;
    let ctd_ixn = "http://ctdbase.org/detail.go?type=relationship&ixnId=";

    let data = fs::read_to_string(&options.input)?;
    let mut ixnset_element = Element::parse(data.as_bytes()).unwrap();

    let model = parse_input(&ixnset_element)?;

    let mut ontology = ontology::set::SetOntology::default();

    for ixn in model.iter() {
        let graph_iri = format!("{}{}", ctd_ixn, ixn.id);
        debug!("graph_iri: {}", graph_iri);
        for (taxon_idx, taxon) in ixn.taxon.iter().enumerate() {
            let ixn_individual = build.named_individual(format!("{}{}#{}", ctd_to_owl_rs::CTDIXN, ixn.id, taxon_idx));
            for (mut codes, actor_individual, mut axioms) in process_ixn(&build, &ixn, &taxon_idx, &chebi_to_mesh_map) {
                if all(&codes, |p| p.as_str() == "w") {
                    axioms.push(Axiom::ObjectPropertyAssertion(ObjectPropertyAssertion::new(
                        build.object_property(format!("{}/{}", OBO, "RO_0002233")).into(),
                        ixn_individual.clone(),
                        actor_individual.clone(),
                    )));
                    axioms.push(Axiom::ClassAssertion(ClassAssertion {
                        ce: ClassExpression::from(build.class(format!("{}/{}", ctd_to_owl_rs::OBO, "CTDI_26"))),
                        i: ixn_individual.clone(),
                    }));
                    add_remnants(&ixn, &build, &mut axioms, &ixn_individual, &taxon);
                }

                // case Interaction("w" :: Nil, _, actors) if actors.forall(_.isInstanceOf[AtomicActor]) =>
                // val (inputs, inputsAxioms) = actors.collect {
                //     case actor: AtomicActor => actor.owl(taxonIndex, meshToCHEBI)
                // }.unzip
                // val cotreatmentAxioms = inputs.map(ixnInd Fact(HasInput, _))
                // Some(ixnInd, inputsAxioms.toSet.flatten ++ cotreatmentAxioms ++ Set(ixnInd Type CoTreatment))

                if all(&codes, |p| p.as_str() == "b") {
                    axioms.push(Axiom::ObjectPropertyAssertion(ObjectPropertyAssertion::new(
                        build.object_property(format!("{}/{}", OBO, "RO_0002233")).into(),
                        ixn_individual.clone(),
                        actor_individual.clone(),
                    )));
                    axioms.push(Axiom::ClassAssertion(ClassAssertion {
                        ce: ClassExpression::from(build.class(format!("{}/{}", ctd_to_owl_rs::OBO, "GO_0005488"))),
                        i: ixn_individual.clone(),
                    }));
                    add_remnants(&ixn, &build, &mut axioms, &ixn_individual, &taxon);
                }

                // case Interaction("b" :: Nil, _, actors) if actors.forall(_.isInstanceOf[AtomicActor]) =>
                // val (inputs, inputsAxioms) = actors.collect {
                //     case actor: AtomicActor => actor.owl(taxonIndex, meshToCHEBI)
                // }.unzip
                // val bindingAxioms = inputs.map(ixnInd Fact(HasInput, _))
                // Some(ixnInd, inputsAxioms.toSet.flatten ++ bindingAxioms ++ Set(ixnInd Type Binding))

                if all(&codes, |p| p.as_str() == "rxn") {

                    //asdf
                }

                // case Interaction("rxn" :: Nil, _, subject :: Interaction(_, innerNode, _) :: Nil) =>
                // val subjectStuff = subject match {
                //     case atomic: AtomicActor =>
                //     val (subjInd, subjAxioms) = atomic.owl(taxonIndex, meshToCHEBI)
                //     val subjProcess = Individual(s"${subjInd.getIRI.toString}-process")
                //     Some((subjProcess, subjAxioms ++ Set(subjProcess Type Process,
                //                                          subjProcess Fact(HasParticipant, subjInd))))
                //     case ixn: Interaction    => interaction(ixn.node, taxonIndex, meshToCHEBI).map(res => (res._1, res._3))
                // }
                // for {
                //     (subjectProcess, subjectAxioms) <- subjectStuff
                //         (_, innerAffector, innerAxioms) <- interaction(innerNode, taxonIndex, meshToCHEBI)
                // } yield {
                //     val axn = (ixnNode \ "axn").head
                //     val relation = processToProcess((axn \ "@degreecode").head.text)
                //     (ixnInd, innerAxioms ++ subjectAxioms ++ Set(
                //         ixnInd Type Process,
                //         innerAffector Fact(PartOf, ixnInd),
                //         subjectProcess Fact(relation, innerAffector)))
                // }

                for axiom in axioms.iter() {
                    ontology.insert(axiom.clone());
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
    prefix_mapping.add_prefix("PMID", ctd_to_owl_rs::PMID).unwrap();
    prefix_mapping.add_prefix("NCBITaxon", ctd_to_owl_rs::NCBITaxon).unwrap();
    prefix_mapping.add_prefix("NCBIGENE", ctd_to_owl_rs::NCBIGENE).unwrap();
    prefix_mapping.add_prefix("MESH", ctd_to_owl_rs::MESH).unwrap();

    let output = fs::File::create(&options.output).ok().unwrap();
    owx::writer::write(&mut io::BufWriter::new(output), &ontology.into(), Some(&prefix_mapping)).ok().unwrap();

    info!("Duration: {}", format_duration(start.elapsed()).to_string());
    Ok(())
}

fn add_remnants(ixn: &IXN, build: &Build, axioms: &mut Vec<Axiom>, ixn_individual: &NamedIndividual, taxon: &Taxon) {
    let pm_ids = ixn.reference.iter().map(|r| format!("{}/{}", ctd_to_owl_rs::PMID, r.pm_id)).collect_vec();

    pm_ids.iter().for_each(|pm_id_iri| {
        axioms.push(Axiom::AnnotationAssertion(AnnotationAssertion::new(
            ixn_individual.clone().into(),
            Annotation { ap: build.annotation_property("http://purl.org/dc/elements/1.1/source"), av: AnnotationValue::Literal(Literal::Simple { literal: pm_id_iri.clone() }) },
        )))
    });
    let organism = build.named_individual(format!("{}#organism", ixn_individual.0));
    axioms.push(Axiom::ClassAssertion(ClassAssertion { ce: ClassExpression::from(build.class(format!("{}{}", ctd_to_owl_rs::NCBITaxon, &taxon.id))), i: organism.clone() }));
    axioms.push(Axiom::AnnotationAssertion(AnnotationAssertion::new(
        organism.clone().into(),
        Annotation {
            ap: build.annotation_property(horned_owl::vocab::RDFS::Label.iri_s()),
            av: AnnotationValue::Literal(Literal::Simple { literal: format!("{}#{}", &taxon.text, &ixn.id) }),
        },
    )));
    axioms.push(Axiom::ObjectPropertyAssertion(ObjectPropertyAssertion::new(
        build.object_property(format!("{}/{}", OBO, "BFO_0000066")).into(),
        ixn_individual.clone(),
        organism.clone(),
    )));
}

fn process_ixn(
    build: &horned_owl::model::Build,
    ixn: &IXN,
    taxon_idx: &usize,
    chebi_to_mesh_map: &collections::HashMap<String, String>,
) -> Vec<(Vec<String>, NamedIndividual, Vec<Axiom>)> {
    let actors = ixn.actors.iter().flat_map(|a| a.flat()).collect_vec();
    // debug!("actors.len(): {:?}", actors.len());
    let mut actions: Vec<AXN> = ixn.axns.clone();
    actions.extend(actors.iter().map(|a| a.axns.clone()).flatten().collect_vec());
    // debug!("actions.len(): {:?}", actions.len());

    debug!("actors: {:?}", actors);
    debug!("actions: {:?}", actions);

    let mapping = actors
        .iter()
        .map(|actor| {
            let codes = actions.iter().filter(|a| a.parent_id == actor.parent_id).map(|a| a.code.to_string()).collect_vec();
            debug!("codes: {:?}", codes);

            let actor_id = actor.id.clone();
            let (actor_class, entity) = match actor.actor_type.as_str() {
                "chemical" => {
                    let chebi_mapping = chebi_to_mesh_map.get(actor_id.as_str());
                    let actor_class = match chebi_mapping {
                        Some(c) => build.class(c.replace("CHEBI:", ctd_to_owl_rs::CHEBI)),
                        None => build.class(actor_id.replace("MESH:", format!("{}", ctd_to_owl_rs::MESH).as_str())),
                    };
                    (actor_class, build.class(format!("{}/{}", ctd_to_owl_rs::OBO, "CHEBI_24431")))
                }
                "gene" => {
                    let actor_class = build.class(actor_id.replace("GENE:", ctd_to_owl_rs::NCBIGENE));
                    (actor_class, build.class(format!("{}/{}", ctd_to_owl_rs::OBO, "SO_0000704")))
                }
                _ => {
                    panic!("should never get here")
                }
            };

            let actor_individual = build.named_individual(format!("{}{}#{}-{}", ctd_to_owl_rs::CTDIXN, actor.parent_id, taxon_idx, actor.position));
            let actor_text = actor.text.as_ref().expect("could not get text");
            let label = format!("{}#{}-{}", actor_text, actor.parent_id, actor.position);

            let mut assertions = Vec::new();

            assertions.push(Axiom::ClassAssertion(ClassAssertion { ce: ClassExpression::from(actor_class.clone()), i: actor_individual.clone() }));
            assertions.push(Axiom::AnnotationAssertion(AnnotationAssertion::new(
                actor_class.clone().0,
                Annotation { ap: build.annotation_property(horned_owl::vocab::RDFS::Label.iri_s()), av: AnnotationValue::Literal(Literal::Simple { literal: actor_text.clone() }) },
            )));
            assertions.push(Axiom::ClassAssertion(ClassAssertion { ce: ClassExpression::from(entity.clone()), i: actor_individual.clone() }));
            assertions.push(Axiom::AnnotationAssertion(AnnotationAssertion::new(
                actor_individual.clone().0,
                Annotation { ap: build.annotation_property(horned_owl::vocab::RDFS::Label.iri_s()), av: AnnotationValue::Literal(Literal::Simple { literal: label.clone() }) },
            )));

            match &actor.form {
                Some(s) => {
                    assertions.push(Axiom::AnnotationAssertion(AnnotationAssertion::new(
                        build.iri(actor_individual.clone()),
                        Annotation { ap: build.annotation_property("http://ctd.example.org/has_form"), av: AnnotationValue::Literal(Literal::Simple { literal: s.to_string() }) },
                    )));
                }
                _ => {}
            };
            (codes, actor_individual, assertions)
        })
        .collect_vec();

    mapping
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
                                debug!("{:?}", ixn);
                            }
                            "reference" => {
                                let reference_pm_id = ixn_child_element.attributes.get("pmid").unwrap();
                                let reference = Reference::new(reference_pm_id.parse::<i32>().unwrap());
                                ixn.reference.push(reference);
                                debug!("{:?}", ixn);
                            }
                            "axn" => {
                                let mut axn = get_axn_from_element(ixn_child_element).unwrap();
                                if ixn.id == axn.parent_id {
                                    ixn.axns.push(axn);
                                }
                                debug!("{:?}", ixn);
                            }
                            "actor" => {
                                let mut actor = get_actor_from_element(ixn_child_element).unwrap();
                                match ixn_child_element.get_text() {
                                    Some(s) => {
                                        actor.text = Some(s.to_string());
                                        if ixn.id == actor.parent_id {
                                            ixn.actors.push(actor);
                                        }
                                        debug!("{:?}", ixn);
                                    }
                                    None => {
                                        if ixn.id == actor.parent_id {
                                            ixn.actors.push(actor);
                                        }
                                        debug!("{:?}", ixn);

                                        for a_node in ixn_child_element.children.iter() {
                                            match a_node {
                                                XMLNode::Element(a_element) => match a_element.name.as_str() {
                                                    "axn" => {
                                                        let mut a_created_axn = get_axn_from_element(a_element).unwrap();
                                                        let found_actor =
                                                            ixn.actors.iter_mut().find(|b_actor| b_actor.id == a_created_axn.parent_id.to_string()).expect("could not find actor");
                                                        found_actor.axns.push(a_created_axn.clone());
                                                        debug!("{:?}", ixn);
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
                                                                debug!("{:?}", ixn);
                                                            }
                                                            None => {
                                                                let found_actor = ixn
                                                                    .actors
                                                                    .iter_mut()
                                                                    .find(|b_actor| b_actor.id == a_created_actor.parent_id.to_string())
                                                                    .expect("could not find actor");
                                                                found_actor.actors.push(a_created_actor.clone());
                                                                debug!("{:?}", ixn);

                                                                for b_node in a_element.children.iter() {
                                                                    match b_node {
                                                                        XMLNode::Element(b_element) => match b_element.name.as_str() {
                                                                            "axn" => {
                                                                                let mut b_created_axn = get_axn_from_element(b_element).unwrap();
                                                                                ixn.actors.iter_mut().for_each(|a_actor| {
                                                                                    let found_actor =
                                                                                        a_actor.actors.iter_mut().find(|b_actor| b_actor.id == b_created_axn.parent_id.to_string());
                                                                                    match found_actor {
                                                                                        Some(a) => a.axns.push(b_created_axn.clone()),
                                                                                        _ => {
                                                                                            debug!("actor not found for: {:?}", b_created_axn);
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
                                                                                                    debug!("actor not found for: {:?}", b_created_actor);
                                                                                                }
                                                                                            }
                                                                                        });
                                                                                    }
                                                                                    None => {}
                                                                                }
                                                                                debug!("{:?}", ixn);
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
    let form_qualifier = element.attributes.get("form_qualifier").cloned();
    let seq_id = element.attributes.get("seq_id").cloned();
    let actor =
        Actor::new(actor_type.to_string(), actor_id.to_string(), actor_position.parse::<i32>().unwrap(), actor_parent_id.parse::<i32>().unwrap(), form, form_qualifier, seq_id);
    Ok(actor)
}

fn get_axn_from_element(element: &Element) -> Result<AXN, Box<dyn error::Error>> {
    let axn_code = element.attributes.get("code").unwrap();
    let axn_degreecode = element.attributes.get("degreecode").unwrap();
    let axn_position = element.attributes.get("position").unwrap();
    let axn_parent_id = element.attributes.get("parentid").unwrap();
    let axn_text = element.get_text().unwrap().to_string();
    let axn = AXN::new(axn_code.into(), axn_degreecode.into(), axn_position.parse::<i32>().unwrap(), axn_parent_id.parse::<i32>().unwrap(), axn_text);
    Ok(axn)
}
