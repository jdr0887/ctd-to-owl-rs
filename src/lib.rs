extern crate env_logger;
extern crate log;
extern crate serde;
extern crate serde_derive;
#[macro_use]
extern crate lazy_static;
extern crate horned_owl;

use std::collections;

pub mod model;

pub const OBO: &str = "http://purl.obolibrary.org/obo";
pub const CTDIXN: &str = "http://ctdbase.org/detail.go?type=relationship&ixnId=";
pub const MESH: &str = "http://id.nlm.nih.gov/mesh/";
pub const CHEBI: &str = "http://purl.obolibrary.org/obo/CHEBI_";
pub const NCBIGENE: &str = "http://identifiers.org/ncbigene/";
pub const PMID: &str = "https://www.ncbi.nlm.nih.gov/pubmed";
pub const NCBI_TAXON: &str = "http://purl.obolibrary.org/obo/NCBITaxon_";
pub const DC: &str = "http://purl.org/dc/elements/1.1";

lazy_static! {
    pub static ref ACTS_UPSTREAM_OF: String = format!("{}/{}", OBO, "RO_0002263");
    pub static ref ACTS_UPSTREAM_OF_POSITIVE_EFFECT: String = format!("{}/{}", OBO, "RO_0004034");
    pub static ref ACTS_UPSTREAM_OF_NEGATIVE_EFFECT: String = format!("{}/{}", OBO, "RO_0004035");
    pub static ref CAUSALLY_UPSTREAM_OF: String = format!("{}/{}", OBO, "RO_0002411");
    pub static ref CAUSALLY_UPSTREAM_OF_POSITIVE_EFFECT: String = format!("{}/{}", OBO, "RO_0002304");
    pub static ref CAUSALLY_UPSTREAM_OF_NEGATIVE_EFFECT: String = format!("{}/{}", OBO, "RO_0002305");
    pub static ref PART_OF: String = format!("{}/{}", OBO, "BFO_0000050");
    pub static ref HAS_PARTICIPANT: String = format!("{}/{}", OBO, "RO_0000057");
    pub static ref HAS_INPUT: String = format!("{}/{}", OBO, "RO_0002233");
    pub static ref INPUT_OF: String = format!("{}/{}", OBO, "RO_0002352");
    pub static ref ENABLES: String = format!("{}/{}", OBO, "RO_0002327");
    pub static ref ENABLED_BY: String = format!("{}/{}", OBO, "RO_0002333");
    pub static ref TRANSPORTS_OR_MAINTAINS_LOCALIZATION_OF: String = format!("{}/{}", OBO, "RO_0002313");
    pub static ref OCCURS_IN: String = format!("{}/{}", OBO, "BFO_0000066");
    pub static ref COTREATMENT: String = format!("{}/{}", OBO, "CTDI_26");
    pub static ref BINDING: String = format!("{}/{}", OBO, "GO_0005488");
    pub static ref PROCESS: String = format!("{}/{}", OBO, "BFO_0000015");
    pub static ref DC_SOURCE: String = format!("{}/{}", DC, "source");
    pub static ref AXN_CODES: Vec<&'static str> = vec![
        "act", "pho", "exp", "myl", "sec", "loc", "clv", "mut", "deg", "spl", "rec", "sta", "met", "oxd", "ubq", "nit", "upt", "red", "alk", "sum", "gyc", "trt", "glc", "csy",
        "upt", "red", "hdx"//, "abu", "ace", "oxd", "fol"
    ];
}

pub fn get_class_map() -> collections::HashMap<String, String> {
    let mut map = collections::HashMap::new();
    map.insert("exp".to_string(), format!("{}/{}", OBO, "GO_0010467"));
    map.insert("w".to_string(), format!("{}/{}", OBO, "CTDI_26"));
    map.insert("rec".to_string(), format!("{}/{}", OBO, "GO_0042221"));
    map.insert("met".to_string(), format!("{}/{}", OBO, "GO_0008152"));
    map.insert("act".to_string(), format!("{}/{}", OBO, "GO_0003674"));
    map.insert("myl".to_string(), format!("{}/{}", OBO, "GO_0032259"));
    map.insert("upt".to_string(), format!("{}/{}", OBO, "CTDI_25"));
    map.insert("imt".to_string(), format!("{}/{}", OBO, "GO_0098657"));
    map.insert("b".to_string(), format!("{}/{}", OBO, "GO_0005488"));
    map.insert("clv".to_string(), format!("{}/{}", OBO, "CTDI_8"));
    map.insert("oxd".to_string(), format!("{}/{}", OBO, "CTDI_20"));
    map.insert("red".to_string(), format!("{}/{}", OBO, "CTDI_21"));
    map.insert("csy".to_string(), format!("{}/{}", OBO, "CTDI_10"));
    map.insert("pho".to_string(), format!("{}/{}", OBO, "GO_0016310"));
    map.insert("loc".to_string(), format!("{}/{}", OBO, "GO_0051179"));
    map.insert("sec".to_string(), format!("{}/{}", OBO, "GO_0046903"));
    map.insert("spl".to_string(), format!("{}/{}", OBO, "GO_0008380"));
    map.insert("ogl".to_string(), format!("{}/{}", OBO, "GO_0006493"));
    map.insert("mut".to_string(), format!("{}/{}", OBO, "CTDI_19"));
    map.insert("trt".to_string(), format!("{}/{}", OBO, "GO_0006810"));
    map.insert("deg".to_string(), format!("{}/{}", OBO, "GO_0009056"));
    map.insert("sta".to_string(), format!("{}/{}", OBO, "CTDI_24"));
    map.insert("ace".to_string(), format!("{}/{}", OBO, "CTDI_2"));
    map.insert("fol".to_string(), format!("{}/{}", OBO, "CTDI_13"));
    map.insert("ubq".to_string(), format!("{}/{}", OBO, "GO_0016567"));
    map.insert("nit".to_string(), format!("{}/{}", OBO, "GO_0017014"));
    map.insert("alk".to_string(), format!("{}/{}", OBO, "CTDI_5"));
    map.insert("sum".to_string(), format!("{}/{}", OBO, "GO_0016925"));
    map.insert("pre".to_string(), format!("{}/{}", OBO, "GO_0018342"));
    map.insert("gyc".to_string(), format!("{}/{}", OBO, "CTDI_15"));
    map.insert("abu".to_string(), format!("{}/{}", OBO, "CTDI_1"));
    map.insert("glc".to_string(), format!("{}/{}", OBO, "GO_0018411"));
    map.insert("hdx".to_string(), format!("{}/{}", OBO, "CTDI_16"));
    map
}
