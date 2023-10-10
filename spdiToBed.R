## SPDI to Y2AVE format converter
## Usage: Rscript spdiToBed.R spdiFile 
##  spdiFile should have only 1 column (canonical SPDI), no header
##  Prints output BED file:
##    chr start end SPDI     (where start = SPDI position, end = SPDI position + length(ref allele), chr = "chr1")
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

## Read the input file
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

spdiBed <- spdiSplit %>%
	transmute(chr=chr, start=position, end=as.numeric(position)+nchar(ReferenceAllele), name=SPDI)

## Write to stdout
write.table(spdiBed, sep='\t', quote=F, col.names=T, row.names=F)
