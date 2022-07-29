library(tidyverse)
library(BSgenome)
library(plyranges)

species_name <- "Crotalus_viridis_viridis"

dfam_lepidosaurs_seq <- readDNAStringSet("seq/compiled_library.fasta")
dfam_lepidosaurs_ranges <- tibble(seqnames = names(dfam_lepidosaurs_seq),
                                  start = 1, end = width(dfam_lepidosaurs_seq)) %>%
  mutate(species = sub(".*@", "", seqnames), short_name = sub(" .*", "", seqnames)) %>%
  as_granges()
base::unique(dfam_lepidosaurs_ranges$species)

dfam_species_ranges <- dfam_lepidosaurs_ranges %>%
  filter(species == species_name) %>%
  mutate(new_name = paste0(species, "_", short_name))

dfam_wo_species_ranges <- dfam_lepidosaurs_ranges[!seqnames(dfam_lepidosaurs_ranges) %in% seqnames(dfam_species_ranges)]

dfam_filtered_seq <- dfam_lepidosaurs_seq[!names(dfam_lepidosaurs_seq) %in% seqnames(dfam_species_ranges)]

dfam_species_seq <- getSeq(dfam_lepidosaurs_seq, dfam_species_ranges)
names(dfam_species_seq) <- dfam_species_ranges$new_name

writeXStringSet(dfam_species_seq, paste0("seq/", species_name, "_solo.fasta"))
writeXStringSet(dfam_filtered_seq, paste0("seq/compiled_library_wo_", species_name, ".fasta"))
