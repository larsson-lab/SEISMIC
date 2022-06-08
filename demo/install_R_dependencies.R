cat('Installing R dependencies\n')
install.packages(c('tidyverse', 'yaml', 'foreach', 'doSNOW', 'stringi', 'inline', 'Rcpp', 'fitdistrplus', 'data.table', 'BiocManager', 'cowplot'))
BiocManager::install("plyranges")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

