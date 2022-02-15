library(dplyr)
library(Biostrings)

#' Detect reads assignment to barcodes for demultiplexing by sample
#' 
#'@reads
#'@barcodeList
#'@max.mismatch int
#'@with.indels boolean
checkBarcodes <- function(reads, barcodesList, max.mismatch=max.mismatch, with.indels=with.indels) {
  
  matrix <- Biostrings::vcountPDict(barcodesList, reads, max.mismatch=max.mismatch, with.indels=with.indels) %>% t()
  matrix <- matrix %>% as.data.frame() %>% dplyr::mutate(across(everything(), ~replace(., . == 1, as.numeric(gsub("V", "", cur_column())))))
  
  if (!any(apply(matrix, 1, function(x) {sum(x != 0)}) > 1)) {
    indexes <- apply(matrix, 1, function(x) {sum(as.numeric(x))})
    indexes[indexes == 0] <- NA
  } else {
    stop("Multiple barcodes match detected. Limit max.mismatch or with.indels arguments.")
  }
  return(indexes)
}