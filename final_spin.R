#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})


# parse input variables
option_list = list(
  make_option(c("-g", "--genome"), type="character", default=NULL,
              help="genome name", metavar="character"),
  make_option(c("-l", "--repeat_library"), type="character", default=NULL,
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

if (is.null(opt$repeat_library)) {
  stop("Repeat library is needed")
} else {
  # alter repeat library name
  opt$repeat_library <- sub(".*/", "", opt$repeat_library)
}

suppressPackageStartupMessages({
  library(BSgenome)
  library(plyranges)
  library(tidyverse)
})

first_washed <- readDNAStringSet(filepath = paste0("data/wash_1_", opt$genome_name))

# make empty holder to add new sequences to
new_consensuses <- DNAStringSet()

to_read <- read_tsv(paste0("data/", opt$genome_name, "_to_rewash_align.txt"), col_names = F, show_col_types = FALSE)

pre_rewash_seq <- readDNAStringSet(paste0("data/needs_rewash_", opt$genome_name))

# read in new sequences and compile
for (i in seq_along(to_read$X1)) {
  if(file.exists(paste0("data/CIAlign/", to_read$X1[i], "_consensus.fasta"))){
    holding_seq <- readDNAStringSet(paste0("data/CIAlign/", to_read$X1[i], "_consensus.fasta"))
    # add old classification
    names(holding_seq) <- names(pre_rewash_seq[sub("#.*", "", names(pre_rewash_seq)) ==
                                                 sub(".fasta", "", sub("rewash_", "", to_read$X1[i]))])
    new_consensuses <- c(new_consensuses, holding_seq)
  } else { # no new consensus use old consensus
    message(paste0("Could not locate data/CIAlign/", to_read$X1[i], "_consensus.fasta"))
    holding_seq <- pre_rewash_seq[sub("#.*", "", names(pre_rewash_seq)) ==
                                    sub(".fasta", "", sub("rewash_", "", to_read$X1[i]))]
    new_consensuses <- c(new_consensuses, holding_seq)
  }
}

# compile first wash with rewashed, write to file
completed_washing <- c(first_washed, new_consensuses)
writeXStringSet(x = completed_washing, filepath = paste0("out/cleaned_library_", opt$repeat_library))
