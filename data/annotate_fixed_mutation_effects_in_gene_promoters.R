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
config$regions_path <- paste0(script_dir, "/gene_promoters_500bp.tsv")
config$mut_effect_path <- paste0(script_dir, "/gene_promoters_500bp_tiled_mut_effect.tsv.gz")
config$fixed_effect <- "na"


str_revcomp <- function(str){
  return(stri_reverse(chartr("TCGA", "AGCT", str)))
}

cat("Reading files\n")
regions.df <- read_tsv(config$regions_path) %>% arrange(gene, start)

regions.gr <- regions.df %>% as_granges()
base_count <- regions.gr %>% as_tibble() %>% group_by(gene) %>% summarise(bases = sum(width))

cat("Tiling regions\n")
tiled_regions.gr <- regions.gr %>%
  tile_ranges(1) %>%
  mutate(gene = rep(base_count$gene, base_count$bases)) %>%
  select(-partition) %>%
  mutate(refnuc = getSeq(BSgenome.Hsapiens.UCSC.hg19, .),
         ref_trinuc = getSeq(BSgenome.Hsapiens.UCSC.hg19, . + 1))

cat("Enumerating bases\n")
tiled_regions.df <- tiled_regions.gr %>% 
  as_tibble() %>% 
  left_join(base_count) %>% 
  group_by(gene) %>% 
  mutate(base_no = ifelse(strand == "+", row_number(), -(row_number() - bases - 1))) %>% 
  mutate(codon_no = 0) %>% # Signifying codons are irrelecant. NA better?
  select(-bases) %>% 
  ungroup()

cat("Making pyr-focused, and annotating effect\n")
pyr <- c("C", "T")
tiled_regions.df <- tiled_regions.df %>%
  mutate(strand = ifelse(refnuc %in% pyr, as.character(strand), chartr("+-", "-+", as.character(strand))),
         ref_trinuc = ifelse(refnuc %in% pyr, ref_trinuc, str_revcomp(ref_trinuc)),
         refnuc = ifelse(refnuc %in% pyr, refnuc, str_revcomp(refnuc))) %>% 
  mutate(C = config$fixed_effect,
         A = config$fixed_effect,
         G = config$fixed_effect,
         T = config$fixed_effect) %>% 
  pivot_longer(c("C","A","G","T"), values_to = "effect", names_to = "mut_to") %>% 
  mutate(effect = ifelse(refnuc == mut_to, "u", effect)) %>% 
  pivot_wider(values_from = effect, names_from = mut_to)


tiled_regions.df %>% select(-width) %>% write_tsv(config$mut_effect_path)
cat("Done!\n")

