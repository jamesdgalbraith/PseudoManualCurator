#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})


# parse input variables
option_list = list(
  make_option(c("-g", "--genome"), type="character", default=NULL,
              help="genome name", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$genome)) {
  stop("Genome name is needed")
} else {
  # set genome names
  opt$genome_name <- sub(".*/", "", opt$genome)
}

suppressPackageStartupMessages({
  library(BSgenome)
  library(plyranges)
  library(tidyverse)
})

# determine new sequences to be read in
to_read <- read_tsv(paste0("data/", opt$genome_name, "_to_align.txt"), col_names = F)

# read in clustered seq
clustered_seq <- readDNAStringSet(paste0("data/clustered_", opt$genome_name))
names(clustered_seq) <- sub(" .*", "", names(clustered_seq))

# make empty holder to add new sequences to
new_consensuses <- DNAStringSet()

# read in new sequences and compile
for (i in seq_along(to_read$X1)) {
  if(file.exists(paste0("data/CIAlign/", to_read$X1[i], "_consensus.fasta"))){
    new_consensuses <- c(new_consensuses, readDNAStringSet(paste0("data/CIAlign/", to_read$X1[i], "_consensus.fasta")))
  }
}

# Rename new sequences to include original classification
names(new_consensuses) <- gsub(".{6}$", "", names(new_consensuses))
compiled_names <- tibble(seqnames = sub("#.*", "", names(clustered_seq)), classified = names(clustered_seq))
better_names <- tibble(seqnames = names(new_consensuses), start = 1, end = width(new_consensuses)) %>%
  inner_join(compiled_names)
names(new_consensuses) <- better_names$classified

# determine those likely fully extended to run through again
to_rerun <-
  inner_join(tibble(seqnames = names(new_consensuses), new_width = width(new_consensuses)),
           tibble(seqnames = names(clustered_seq), old_width = width(clustered_seq))) %>%
  filter(new_width >= 1400 + old_width)

# compile completed (not to rerun) and not improved seq
not_improved_seq <- clustered_seq[!names(clustered_seq) %in% names(new_consensuses)]
completed_seq <- new_consensuses[!names(new_consensuses) %in% to_rerun$seqnames]
writeXStringSet(x = c(not_improved_seq, completed_seq), filepath = paste0("data/wash_1_", opt$genome_name))

# get sequences to rewash
to_rerun_seq <- new_consensuses[names(new_consensuses) %in% to_rerun$seqnames]
writeXStringSet(to_rerun_seq, paste0("data/needs_rewash_", opt$genome_name))

