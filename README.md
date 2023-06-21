# Demo of using SPDI tools for IGVF variant lists

Here I provide a demo for some of the NCBI tools for adding SPDI
unique IDs to variant lists. Briefly the benefits of SPDI:

* resolves indel ambiguity
* validates, e.g. incorrect reference allele specification
* human readable, unique ID with broad usage across consortia
  (dbSNP, ClinVar, etc.), NCBI support including API and toolkits

See the following for a tutorial on NCBI tools:

<https://github.com/ncbi/dbsnp/tree/master/tutorials/Variation%20Services>

Here I (Mike) have added an option `SPDI` to the NCBI scripts that
will convert from...

```
chrom : position (1-based) : ref : alt
```

(remove spaces)

...to a unique SPDI (0-based and validated).

## IGVF use case:

Suppose we have a variant list that looks like "1_25253604_hg38_G_A"
(this is customizable with simple arguments within `make_spdi_list.R`).

```
> Rscript make_spdi_list.R variants.txt
> head -100 spdi_for_batch_processing.txt > spdi_100.txt
> python spdi_batch.py -i spdi_100.txt -t SPDI

NC_000001.11:25253603:G:A	NC_000001.11:25253603:G:A
NC_000001.11:25336579:C:G	NC_000001.11:25336579:C:G
NC_000001.11:25341419:G:A	NC_000001.11:25341419:G:A
NC_000001.11:25341834:C:T	NC_000001.11:25341834:C:T
NC_000001.11:25342222:T:C	NC_000001.11:25342222:T:C
NC_000001.11:25348293:C:T	NC_000001.11:25348293:C:T
...
```

## Example with HGVS

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
