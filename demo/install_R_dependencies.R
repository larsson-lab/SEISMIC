cat('Installing R dependencies\n')
install.packages('tidyverse', 'yaml', 'foreach', 'doParallel', 'stringi', 'inline', 'Rcpp', 'fitdistrplus', 'data.table', 'BiocManager')
BiocManager::install("plyranges")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
