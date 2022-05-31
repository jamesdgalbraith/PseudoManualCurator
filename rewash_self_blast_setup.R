#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})


# parse input variables
option_list = list(
  make_option(c("-g", "--genome"), type="character", default=NULL,
              help="genome name", metavar="character"),
  make_option(c("-l", "--library"), type="character", default=NULL,
              help="genome name", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$genome)) {
  stop("Path to genome is needed")
} else {
  # set genome names
  opt$genome_name <- sub(".*/", "", opt$genome)
}

suppressPackageStartupMessages({
  library(BSgenome)
  library(plyranges)
  library(tidyverse)
})

flank5 = 1500
flank3 = 1500

blast_out <- read_tsv(file = paste0("data/", opt$genome_name, "_rewash_search.out"),
                      col_names = c("qseqid", "seqnames", "pident", "length", "qstart", "qend",
                                    "qlen", "sstart", "send", "slen", "evalue", "bitscore"))

blast_out_trimmed <- blast_out %>%
  dplyr::filter(
    pident >= 80,
    length >= 0.5 * qlen
  ) %>%
  dplyr::group_by(qseqid) %>%
  dplyr::arrange(-bitscore) %>%
  dplyr::slice(1:20) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(strand = base::ifelse(sstart < send, "+", "-"),
                start = base::ifelse(sstart < send, sstart - flank5, send - flank3),
                end = base::ifelse(sstart > send, sstart + flank3, send + flank5),
                start = ifelse(start < 1, 1, start),
                end = ifelse(end > slen, slen, end),
                seqnames = paste0(qseqid, "___", seqnames)) %>%
  dplyr::select(seqnames, qseqid, start, end, strand) %>%
  plyranges::as_granges()

blast_out_merged <- as_tibble(plyranges::reduce_ranges(blast_out_trimmed)) %>%
  mutate(seqnames = as.character(seqnames)) %>%
  separate(seqnames, into = c("qseqid", "seqnames"), sep = "___") %>%
  dplyr::select(-strand, -width) %>%
  plyranges::as_granges()

blast_out_tbl <- tibble(qseqid = names(table(blast_out_merged$qseqid)), n = as.integer(table(blast_out_merged$qseqid))) %>%
  filter(n > 2)

blast_out_ignored_tbl <- tibble(qseqid = names(table(blast_out_merged$qseqid)), n = as.integer(table(blast_out_merged$qseqid))) %>%
  filter(n <= 2)

genome_seq <- Biostrings::readDNAStringSet(filepath = opt$genome)
names(genome_seq) <- sub(" .*", "", names(genome_seq))

consensus_seq <- Biostrings::readDNAStringSet(filepath = opt$library)
names(consensus_seq) <- sub(" .*", "", names(consensus_seq))

for(i in seq_along(blast_out_tbl$qseqid)){
  
  align_ranges <- blast_out_merged[blast_out_merged$qseqid == blast_out_tbl$qseqid[i]]
  align_ranges <- GenomicRanges::reduce(align_ranges, ignore.strand=T)
  align_seq <- BSgenome::getSeq(genome_seq, align_ranges)
  names(align_seq) <- base::paste0(GenomeInfoDb::seqnames(align_ranges), ":", IRanges::ranges(align_ranges))
  
  align_seq <- c(consensus_seq[names(consensus_seq) == blast_out_tbl$qseqid[i]], align_seq)
  
  Biostrings::writeXStringSet(x = align_seq, filepath = paste0("data/initial_seq/rewash_", opt$genome_name, "_",
                                                               sub("#.*", "", blast_out_tbl$qseqid[i]), ".fasta"))
  
}

