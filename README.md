# Demo of using SPDI tools for IGVF variant lists

This repository contains some tools/scripts for adding SPDI
unique IDs to variant lists. Briefly the benefits of SPDI:

* resolves indel ambiguity
* validates, e.g. incorrect reference allele specification
* human readable, unique ID with broad usage across consortia
  (dbSNP, ClinVar, etc.), NCBI support including API and toolkits

There is a script `spdi_batch.py` that leverages the NCBI API to provide 
SPDI. There is also an R script `jon_rosen_make_spdi.R` that provides
a faster version that doesn't make use of an API but normalizes alleles
locally (and has been checked across 1000s of indels). We provide the latter
although for ultimate certainty for aligning with NCBI you should check 
with the former. 

The scripts here will provide the **Canonical SPDI** for variants. This 
means that it is the *unique representation*, particularly important for 
indels which can have arbitrary allelic representation. You must use the 
Canonical SPDI to refer to indels in order to allow comparison of variants 
across centers. The Canonical SPDI will appear in the right column after
running the script. It will not give a warning or error, but will fix the
allelic representation if it is not canonical, so you should check to see
if the left and right columns differ for any rows.

See the following for a tutorial on NCBI tools:

<https://github.com/ncbi/dbsnp/tree/master/tutorials/Variation%20Services>

## Additions:

* I have added an example script `make_spdi_list.R` for converting from
arbitrary variant specification (1-based) to SPDI input (0-based) for
the NCBI Variant Services.

* I have added an option `SPDI` to the NCBI's original script
  `spdi_batch.py` that will convert from 0-based input:

```
chrom:position:ref:alt
```

...to a unique SPDI (0-based and validated).

## IGVF use case:

Suppose we have a variant list that looks like "1_25253604_hg38_G_A"
(1-based positions, separated by underscore), in the file
`variants.txt`.

Note that this is easily customizable with arguments within
`make_spdi_list.R`.

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

**What does `warnings` mean?** This typically means that you have
mis-specified the reference allele of hg38.

You can check here (don't forget that the above positions are 0-based
while the genome browser is 1-based):

<https://www.ncbi.nlm.nih.gov/genome/gdv/browser/genome/?id=GCF_000001405.40>

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
