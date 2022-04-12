library(plyranges)
library(tidyverse)
library(stringi)
library(BSgenome.Hsapiens.UCSC.hg19)

thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
script_dir <- dirname(thisFile())

config <- list()
config$single_transcript_gene_list_path <- paste0(script_dir, "/gene_cds.tsv")
config$codon_table_path <- paste0(script_dir, "/codon_table.tsv")
config$transcript_mut_effect_path <- paste0(script_dir, "/gene_cds_tiled_mut_effect.tsv.gz")


cat("Reading files\n")
codon_table.df <- read_tsv(config$codon_table_path)
exons.df <- read_tsv(config$single_transcript_gene_list_path) %>% arrange(gene, start)

# Temporary limit for testing purposes
#exons.df <- exons.df %>% filter(gene %in% c("A1BG", "A1CF", "A2M", "A2ML1", "A3GALT2", "A4GALT"))


exons.gr <- exons.df %>% as_granges()
codon_count <- exons.gr %>% as_tibble() %>% group_by(gene) %>% summarise(bases = sum(width)) %>% mutate(codons = bases / 3)

cat("Tiling exons\n")
tiled_exons.gr <- exons.gr %>%
  tile_ranges(1) %>%
  mutate(gene = rep(codon_count$gene, codon_count$bases)) %>%
  select(-partition) %>%
  mutate(refnuc = getSeq(BSgenome.Hsapiens.UCSC.hg19, .),
         ref_trinuc = getSeq(BSgenome.Hsapiens.UCSC.hg19, . + 1))

cat("Enumerating codons and bases\n")
tiled_exons.df <- tiled_exons.gr %>% 
  as_tibble() %>% 
  left_join(codon_count %>% select(-codons)) %>% 
  group_by(gene) %>% 
  mutate(base_no = ifelse(strand == "+", row_number(), -(row_number() - bases - 1))) %>% 
  mutate(codon_no = ceiling(base_no / 3)) %>% 
  select(-bases)

cat("Getting sequences\n")
gene_seqs.df <- tiled_exons.df %>%
  group_by(gene) %>%
  summarise(seq = ifelse(unique(strand) == "+", paste0(refnuc, collapse = ""), reverse(paste0(refnuc, collapse = ""))))

cat("Finding mutation effects for each base\n")
gene_codons.df <-
  gene_seqs.df %>%
  mutate(., codon = lapply(1:nrow(.), function(x) strsplit(gene_seqs.df$seq[x], "(?<=.{3})", perl = TRUE)[[1]])) %>%
  select(-seq) %>%
  unnest(cols = c(codon)) %>% 
  group_by(gene) %>% 
  mutate(codon_no = row_number()) %>% 
  ungroup() %>% 
  dplyr::slice(rep(1:n(), each = 3)) %>% 
  mutate(codon_base_no = rep(1:3, n()/3)) %>% 
  dplyr::slice(rep(1:n(), each = 4)) %>% 
  mutate(mut_to = rep(c("A", "C", "G", "T"), n()/4)) %>% 
  mutate(codon1 = substr(codon, 1, codon_base_no - 1),
         codon2 = substr(codon, codon_base_no + 1, 3)) %>% 
  mutate(mut_codon = paste0(codon1, mut_to, codon2)) %>% 
  select(-codon1, -codon2) %>% 
  left_join(codon_table.df) %>% 
  left_join(codon_table.df %>% dplyr::rename(mut_codon = codon, mut_aa = aa)) %>% 
  mutate(effect = case_when(codon == mut_codon ~ "u", aa == mut_aa ~ "s", mut_aa == "STOP" ~ "n", TRUE ~ "m"))
# u - unmutated
# s - synonymous
# n - nonsense
# m - missense

cat("Sorting mutation effects and making pyrimidine-focused\n")
str_revcomp <- function(str){
  return(stri_reverse(chartr("TCGA", "AGCT", str)))
}
pyr <- c("C", "T")
gr_sorted_effects.df <-
  gene_codons.df %>%
  mutate(refnuc = substr(codon, codon_base_no, codon_base_no)) %>%
  mutate(strand_switch = ifelse(refnuc %in% pyr, FALSE, TRUE)) %>% 
  mutate(mut_to = ifelse(refnuc %in% pyr, mut_to, str_revcomp(mut_to)),
         #codon = ifelse(refnuc %in% pyr, codon, str_revcomp(codon)), # Don't need the codon, and it's confusing after changing to pyr-focused
         refnuc = ifelse(refnuc %in% pyr, refnuc, str_revcomp(refnuc))) %>% 
  select(gene, codon_no, codon_base_no, mut_to, effect, strand_switch) %>% 
  pivot_wider(names_from = mut_to, values_from = effect) %>% 
  left_join(exons.df %>% select(strand, gene) %>% unique()) %>% 
  left_join(codon_count %>% select(-bases)) %>% 
  mutate(plus_strand_codon_order = ifelse(strand == "+", codon_no, -(codon_no - codons - 1)),
         plus_strand_codon_base_order = ifelse(strand == "+", codon_base_no, -(codon_base_no - 4))) %>% 
  arrange(gene, plus_strand_codon_order, plus_strand_codon_base_order)

cat("Adding mutation effects to tiled exons\n")
tiled_exons.df <- tiled_exons.df %>%
  ungroup() %>%
  bind_cols(gr_sorted_effects.df %>% select(C,A,G,T, strand_switch)) %>%
  mutate(strand = ifelse(strand_switch, chartr("+-", "-+", as.character(strand)), as.character(strand)),
         refnuc = ifelse(strand_switch, str_revcomp(refnuc), refnuc),
         ref_trinuc = ifelse(strand_switch, str_revcomp(ref_trinuc), ref_trinuc))

tiled_exons.df %>% select(-width, -strand_switch) %>% write_tsv(config$transcript_mut_effect_path)
cat("Done!\n")

