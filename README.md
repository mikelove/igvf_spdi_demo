# Demo of using SPDI tools for IGVF variant lists

Here I provide a demo for some of the NCBI tools for adding SPDI
unique IDs to variant lists. Briefly the benefits of SPDI:

* resolves indel ambiguity
* validates unique IDs, e.g. incorrect alt/ref
* human readable, unique ID with broad usage across consortia
  (dbSNP, ClinVar, etc.), NCBI support including API and toolkits

See the following for a tutorial on NCBI tools:

<https://github.com/ncbi/dbsnp/tree/master/tutorials/Variation%20Services>

## Example:

```
> python spdi_batch.py -i test.txt -t HGVS
NC_000021.9:g.25716261G>A	NC_000021.9:25716260:G:A
ERROR: status code = 400
NC_000021.9:g.25716536_25716537insAT	NC_000021.9:25716536::AT
NC_000021.9:g.25716557del	NC_000021.9:25716556:TTTT:TTT
NC_000021.9:g.25716558del	NC_000021.9:25716556:TTTT:TTT
NC_000021.9:g.25716559del	NC_000021.9:25716556:TTTT:TTT
NC_000021.9:g.25716560del	NC_000021.9:25716556:TTTT:TTT
NC_000021.9:g.25716557dup	NC_000021.9:25716556:TTTT:TTTTT
NC_000021.9:g.25716558dup	NC_000021.9:25716556:TTTT:TTTTT
NC_000021.9:g.25716559dup	NC_000021.9:25716556:TTTT:TTTTT
NC_000021.9:g.25716560dup	NC_000021.9:25716556:TTTT:TTTTT
```
