## SPDI to Y2AVE format converter
## Usage: Rscript spdiToY2AVE.R spdiFile 
##  spdiFile should have only 1 column (canonical SPDI), no header
##  Prints output table containing Y2AVE columns:
##    chrRefSeqID chr position ReferenceAllele AlternativeAllele SPDI
## Jesse Engreitz - Oct 9, 2023

suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  cat("Script requires one argument: Input SPDI file\n")
  quit("no")
}
infile <- args[1]

# no scientific notation
options(scipen=999)

## Read input file
spdi <- read.table(infile, stringsAsFactors=F)
colnames(spdi) <- "SPDI"

# Map between GRCh38 molecule name and RefSeq sequence name
RefSeq_start <- rep("NC_0000", 23)
RefSeq_mid <- c(paste0("0", as.character(1:9)),
  as.character(10:23))
RefSeq_end <- as.character(c(11, 12, 12, 12, 10, 12, 14, 11, 12, 11,
  10, 12, 11, 9, 10, 10, 11, 10, 10, 11, 9, 11, 11))
RefSeq <- paste0(RefSeq_start, RefSeq_mid, ".", RefSeq_end) 
RefSeq_map <- data.frame(chr = paste0("chr", c(1:22, "X")), 
                NC = RefSeq)

## Process and format
spdiSplit <- as.data.frame(do.call(rbind, strsplit(spdi$SPDI, ":")), stringsAsFactors=F)
colnames(spdiSplit) <- c("chrRefSeqID","position","ReferenceAllele","AlternativeAllele")
spdiSplit$SPDI <- spdi$SPDI
spdiSplit <- merge(spdiSplit, RefSeq_map, by.x="chrRefSeqID", by.y="NC", all.x=TRUE)
stopifnot(all(!is.na(spdiSplit$chr)))
if (any(spdiSplit$ReferenceAllele == ""))
  spdiSplit$ReferenceAllele[spdiSplit$ReferenceAllele == ""] <- "-"
if (any(spdiSplit$AlternativeAllele == ""))
  spdiSplit$AlternativeAllele[spdiSplit$AlternativeAllele == ""] <- "-"

spdiSplit <- spdiSplit[,c("chrRefSeqID","chr","position","ReferenceAllele","AlternativeAllele","SPDI")]

## Write to stdout
write.table(spdiSplit, sep='\t', quote=F, col.names=T, row.names=F)
