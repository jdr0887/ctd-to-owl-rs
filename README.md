# ctd-to-owl-rs

This project is an attempt at using Rust to speed up XML parsing & processing of [CTD](http://ctdbase.org) data to generate an OWL ontology.  The original code this project is based on was written in Scala.

The [ctd-to-owl](https://github.com/balhoff/ctd-to-owl) project runs in just under 50 minutes and requires >120GB of memory, while this Rust based implementation runs in minutes (2:04 is the best I have witnessed) while using <30GB of memory.

Example Build/Usage:
```
ctd-to-owl-rs$ cargo build --release
ctd-to-owl-rs$ RUST_LOG=info ./target/release/ctd-to-owl -i <some_dir>/CTD_chem_gene_ixns_structured.xml -o <some_dir> -c <some_dir>/chebi_mesh.tsv
```
