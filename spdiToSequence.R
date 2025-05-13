# Script to take in SPDI and output allele sequences (with flanking bp)
# Michael Love
# May 13 2025

# Rscript spdiToSequence.R test_data/spdi_100.txt 10 test_out.txt

# arguments:
# 1 - file with list of SPDIs
# 2 - number of basepairs to flank on either side
# 3 - output filename (3 cols: spdi, ref, alt)
args <- commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(plyranges)
  library(BSgenome.Hsapiens.NCBI.GRCh38)
})

spdi <- scan(args[1], what="char")
dat <- tibble(spdi)

gr <- dat |>
  separate(spdi, into=c("seqnames","start","ref","alt"), sep=":") |>
  mutate(seqnames = str_remove(seqnames, "NC_(0+)"),
         seqnames = str_remove(seqnames, "\\..*"),
         seqnames = case_when(
           seqnames == "23" ~ "X",
           seqnames == "24" ~ "Y",
           seqnames == "12920" ~ "MT",
           TRUE ~ seqnames
         ),
         start = as.numeric(start) + 1,
         width = nchar(ref)) |>
  as_granges()

ref_seq <- getSeq(Hsapiens, gr)
stopifnot(all.equal(gr$ref, as.character(ref_seq)))

add_to_side <- as.numeric(args[2]) # amount to add to each side

left_seq <- gr |>
  flank_left(add_to_side) %>%
  getSeq(Hsapiens, .)

right_seq <- gr |>
  flank_right(add_to_side) %>%
  getSeq(Hsapiens, .)

ref <- paste(left_seq, gr$ref, right_seq, sep="")
alt <- paste(left_seq, gr$alt, right_seq, sep="")

write.table(cbind(spdi, ref, alt), file=args[3], quote=FALSE, row.names=FALSE, col.names=FALSE)
message(paste("Wrote",length(spdi),"items"))
