extern crate serde;
#[macro_use]
extern crate serde_derive;
#[macro_use]
extern crate log;
extern crate env_logger;
#[macro_use]
extern crate lazy_static;
extern crate horned_owl;

use std::collections;
use std::error;

pub mod model;

pub const OBO: &str = "http://purl.obolibrary.org/obo";
pub const CTDIXN: &str = "http://ctdbase.org/detail.go?type=relationship&ixnId=";
pub const MESH: &str = "http://id.nlm.nih.gov/mesh/";
pub const CHEBI: &str = "http://purl.obolibrary.org/obo/CHEBI_";
pub const NCBIGENE: &str = "http://identifiers.org/ncbigene/";
pub const PMID: &str = "https://www.ncbi.nlm.nih.gov/pubmed";
pub const NCBITaxon: &str = "http://purl.obolibrary.org/obo/NCBITaxon_";

pub fn get_class_map(build: &horned_owl::model::Build) -> Result<collections::HashMap<String, horned_owl::model::Class>, Box<dyn error::Error>> {
    let mut map = collections::HashMap::new();
    map.insert("exp".to_string(), build.class(format!("{}/{}", OBO, "GO_0010467")));
    map.insert("w".to_string(), build.class(format!("{}/{}", OBO, "CTDI_26")));
    map.insert("rec".to_string(), build.class(format!("{}/{}", OBO, "GO_0042221")));
    map.insert("met".to_string(), build.class(format!("{}/{}", OBO, "GO_0008152")));
    map.insert("act".to_string(), build.class(format!("{}/{}", OBO, "GO_0003674")));
    map.insert("myl".to_string(), build.class(format!("{}/{}", OBO, "GO_0032259")));
    map.insert("upt".to_string(), build.class(format!("{}/{}", OBO, "CTDI_25")));
    map.insert("imt".to_string(), build.class(format!("{}/{}", OBO, "GO_0098657")));
    map.insert("b".to_string(), build.class(format!("{}/{}", OBO, "GO_0005488")));
    map.insert("clv".to_string(), build.class(format!("{}/{}", OBO, "CTDI_8")));
    map.insert("oxd".to_string(), build.class(format!("{}/{}", OBO, "CTDI_20")));
    map.insert("red".to_string(), build.class(format!("{}/{}", OBO, "CTDI_21")));
    map.insert("csy".to_string(), build.class(format!("{}/{}", OBO, "CTDI_10")));
    map.insert("pho".to_string(), build.class(format!("{}/{}", OBO, "GO_0016310")));
    map.insert("loc".to_string(), build.class(format!("{}/{}", OBO, "GO_0051179")));
    map.insert("sec".to_string(), build.class(format!("{}/{}", OBO, "GO_0046903")));
    map.insert("spl".to_string(), build.class(format!("{}/{}", OBO, "GO_0008380")));
    map.insert("ogl".to_string(), build.class(format!("{}/{}", OBO, "GO_0006493")));
    map.insert("mut".to_string(), build.class(format!("{}/{}", OBO, "CTDI_19")));
    map.insert("trt".to_string(), build.class(format!("{}/{}", OBO, "GO_0006810")));
    map.insert("deg".to_string(), build.class(format!("{}/{}", OBO, "GO_0009056")));
    map.insert("sta".to_string(), build.class(format!("{}/{}", OBO, "CTDI_24")));
    map.insert("ace".to_string(), build.class(format!("{}/{}", OBO, "CTDI_2")));
    map.insert("fol".to_string(), build.class(format!("{}/{}", OBO, "CTDI_13")));
    map.insert("ubq".to_string(), build.class(format!("{}/{}", OBO, "GO_0016567")));
    map.insert("nit".to_string(), build.class(format!("{}/{}", OBO, "GO_0017014")));
    map.insert("alk".to_string(), build.class(format!("{}/{}", OBO, "CTDI_5")));
    map.insert("sum".to_string(), build.class(format!("{}/{}", OBO, "GO_0016925")));
    map.insert("pre".to_string(), build.class(format!("{}/{}", OBO, "GO_0018342")));
    map.insert("gyc".to_string(), build.class(format!("{}/{}", OBO, "CTDI_15")));
    map.insert("abu".to_string(), build.class(format!("{}/{}", OBO, "CTDI_1")));
    map.insert("glc".to_string(), build.class(format!("{}/{}", OBO, "GO_0018411")));
    map.insert("hdx".to_string(), build.class(format!("{}/{}", OBO, "CTDI_16")));
    Ok(map)
}
