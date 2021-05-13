use itertools::Itertools;
use serde::Deserialize;

#[derive(Deserialize, Debug, PartialEq, Clone)]
pub struct Taxon {
    pub id: i32,
    pub text: String,
}

impl Taxon {
    pub fn new(id: i32, text: String) -> Taxon {
        Taxon { id, text }
    }
}

#[derive(Deserialize, Debug, PartialEq, Clone)]
pub struct Reference {
    #[serde(rename(deserialize = "pmid"))]
    pub pm_id: i32,
}

impl Reference {
    pub fn new(pm_id: i32) -> Reference {
        Reference { pm_id }
    }
}

#[derive(Deserialize, Debug, PartialEq, Clone)]
pub struct AXN {
    pub code: String,
    #[serde(rename(deserialize = "degreecode"))]
    pub degree_code: char,
    pub position: i8,
    #[serde(rename(deserialize = "parentid"))]
    pub parent_id: i32,
    pub text: String,
}

impl AXN {
    pub fn new(code: String, degree_code: char, position: i8, parent_id: i32, text: String) -> AXN {
        AXN { code, degree_code, position, parent_id, text }
    }
}

#[derive(Deserialize, Debug, PartialEq, Clone)]
pub struct Actor {
    #[serde(rename(deserialize = "type"))]
    pub actor_type: String,
    pub id: String,
    pub position: i8,
    #[serde(rename(deserialize = "parentid"))]
    pub parent_id: i32,
    pub form: Option<String>,
    pub form_qualifier: Option<String>,
    pub seq_id: Option<String>,
    pub text: Option<String>,
    #[serde(rename(deserialize = "axn"))]
    pub axns: Vec<AXN>,
    pub actors: Vec<Actor>,
    // #[serde(rename(deserialize = "actor"))]
    // pub nested_actors: Option<Vec<Actor>>,
}

impl Actor {
    pub fn new(actor_type: String, id: String, position: i8, parent_id: i32, form: Option<String>, form_qualifier: Option<String>, seq_id: Option<String>) -> Actor {
        Actor { actor_type, id, position, parent_id, form, form_qualifier, seq_id, text: None, axns: Vec::new(), actors: Vec::new() }
    }
    pub fn flat(&self) -> Vec<&Actor> {
        std::iter::once(self).chain(self.actors.iter().flat_map(|c| c.flat())).collect_vec()
    }
}

#[derive(Deserialize, Debug, PartialEq, Clone)]
pub struct IXN {
    pub id: i32,
    pub taxon: Vec<Taxon>,
    pub reference: Vec<Reference>,
    pub axns: Vec<AXN>,
    pub actors: Vec<Actor>,
}

impl IXN {
    pub fn new() -> IXN {
        IXN { id: 0, taxon: Vec::new(), reference: Vec::new(), axns: Vec::new(), actors: Vec::new() }
    }
}

#[derive(Deserialize, Debug, PartialEq, Clone)]
pub struct IXNSet {
    #[serde(rename(deserialize = "ixn"))]
    pub ixns: Vec<IXN>,
}

impl IXNSet {
    pub fn new() -> IXNSet {
        IXNSet { ixns: Vec::new() }
    }
    pub fn ixn_mut(&mut self) -> &mut Vec<IXN> {
        &mut self.ixns
    }
}

pub struct Interaction {
    pub codes: Vec<String>,
    pub actors: Vec<Actor>,
}

impl Interaction {
    pub fn new(codes: Vec<String>, actors: Vec<Actor>) -> Interaction {
        Interaction { codes, actors }
    }
}
