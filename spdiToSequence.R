library(dplyr)
library(tidyr)
library(stringr)
library(plyranges)

spdi <- scan("test_data/spdi_100.txt", what="char")
dat <- tibble(spdi)

gr <- dat |>
  separate(spdi, into=c("seqnames","start","ref","alt"), sep=":") |>
  mutate(seqnames = str_remove(seqnames, "NC_(0+)"),
         seqnames = str_remove(seqnames, "\\..*"),
         start = as.numeric(start) + 1,
         width = nchar(ref)) |>
  as_granges()

suppressPackageStartupMessages(library(BSgenome.Hsapiens.NCBI.GRCh38))
ref_seq <- getSeq(Hsapiens, gr)
stopifnot(all.equal(gr$ref, as.character(ref_seq)))

add_to_side <- 5 # amount to add to each side

left_seq <- gr |>
  flank_left(add_to_side) %>%
  getSeq(Hsapiens, .)

right_seq <- gr |>
  flank_right(add_to_side) %>%
  getSeq(Hsapiens, .)

