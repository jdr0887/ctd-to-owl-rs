# ctd-to-owl-rs

This project is an attempt at using Rust to speed up XML parsing & processing of [CTD](http://ctdbase.org) data written by Jim Balhoff using Scala.

The original [ctd-to-owl](https://github.com/balhoff/ctd-to-owl) project runs in just under 50 minutes and requires >120MB of memory, while this Rust based implementation runs in just over 8 minutes while using <30MB of memory.

