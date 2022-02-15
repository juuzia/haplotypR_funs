#' deplexBySample function with extension to allow mismatches and indels
#' and additionally search barcodes of both orientations in each FASTQ file
#' 
#'@fastqFileFwd
#'@fastqFileRev
#'@barcodeFileFws
#'@barcodeFileRev
#'@outputDir
#'@adapterFwd
#'@adapterRev
#'@max.mismatch
#'@with.indels
deplexSampleExtended <- function (fastqFileFwd, fastqFileRev, barcodeFileFwd, barcodeFileRev, 
                    outputDir, adapterFwd = NULL, adapterRev = NULL, max.mismatch = 0, 
                    with.indels = F, progressReport = message) 
{
  require(Biostrings)
  require(ShortRead)
  barcodesFwd <- Biostrings::readDNAStringSet(barcodeFileFwd)
  barcodesFwdLength <- unique(width(barcodesFwd))
  if (length(barcodesFwdLength) > 1) 
    stop("Barcodes length must have equal length.")
  barcodesRev <- Biostrings::readDNAStringSet(barcodeFileRev)
  barcodesRevLength <- unique(width(barcodesRev))
  if (length(barcodesRevLength) > 1) 
    stop("Barcodes length must have equal length.")
  of <- list.files(path = outputDir, pattern = names(barcodesFwd))
  if (length(of) > 0) 
    stop("Output directory must be empty. Found existing files: ", 
         paste(of, collapse = ", "))
  if (!is.function(progressReport)) 
    progressReport <- message
  f1 <- FastqStreamer(fastqFileFwd)
  f2 <- FastqStreamer(fastqFileRev)
  mode <- "w"
  countReads <- 0
  sumDemultiplex <- character(0)
  msg <- paste("Processing file", basename(fastqFileFwd), 
               "and", basename(fastqFileRev), ":")
  progressReport(detail = msg)
  while (length(sr1 <- yield(f1)) > 0 & length(sr2 <- yield(f2)) > 
         0) {
    if (length(sr1) != length(sr2)) {
      warning("Fastq files have unequal lengths, it is likely one is incomplete")
      if (length(sr1) < length(sr2)) 
        sr2 <- sr2[seq_along(sr1)]
      else sr1 <- sr1[seq_along(sr2)]
    }
    countReads <- countReads + length(sr1)
    sr1_trim <- narrow(sr1, start = 1, width = barcodesFwdLength)
    sr2_trim <- narrow(sr2, start = 1, width = barcodesRevLength)
    
    ### Fwd barcode in R1 & Rev barcode in R2 ###
    # R1
    idxFwd_sr1 <- checkBarcodes(reads = sread(sr1_trim),
                                barcodesList = barcodesFwd,
                                max.mismatch = max.mismatch,
                                with.indels = with.indels)
    BFwd_sr1 <- names(barcodesFwd)[as.numeric(idxFwd_sr1)]
    
    # R2
    idxRev_sr2 <- checkBarcodes(reads = sread(sr2_trim),
                                barcodesList = barcodesRev,
                                max.mismatch = max.mismatch,
                                with.indels = with.indels)
    BRev_sr2 <- names(barcodesRev)[as.numeric(idxRev_sr2)]

    # Write sequences to FASTQ without undetected reads
    idx <- paste(BFwd_sr1, BRev_sr2, sep = "-")
    
    sumDemultiplex <- c(sumDemultiplex, idx[idx != "NA-NA"])
    sr1_lst <- split(sr1, idx)
    sr1_lst <- sr1_lst[names(sr1_lst) != "NA-NA"]
    sr2_lst <- split(sr2, idx)
    sr2_lst <- sr2_lst[names(sr2_lst) != "NA-NA"]
    
    lapply(seq_along(sr1_lst), function(i) {
      outFile <- file.path(outputDir, names(sr1_lst[i]))
      writeFastq(sr1_lst[[i]], file = paste(outFile, "_R1.fastq.gz", 
                                            sep = ""), mode = mode, compress = T)
      writeFastq(sr2_lst[[i]], file = paste(outFile, "_R2.fastq.gz", 
                                            sep = ""), mode = mode, compress = T)
    })
    
    ### Fwd barcode in R2 & Rev barcode in R1 ###
    mode <- "a"
    # R2
    idxFwd_sr2 <- checkBarcodes(reads = sread(sr2_trim),
                                barcodesList = barcodesFwd,
                                max.mismatch = max.mismatch,
                                with.indels = with.indels)
    BFwd_sr2 <- names(barcodesFwd)[as.numeric(idxFwd_sr2)]
    
    # R1
    idxRev_sr1 <- checkBarcodes(reads=sread(sr1_trim), barcodesList=barcodesRev, max.mismatch=max.mismatch, with.indels=with.indels)
    BRev_sr1 <- names(barcodesRev)[as.numeric(idxRev_sr1)]

    # Append sequences to FASTQ without undetected reads
    idxOdd <- paste(BFwd_sr2, BRev_sr1, sep = "-")
    
    sumDemultiplex <- c(sumDemultiplex, idxOdd[idxOdd != "NA-NA"])
    sr1_lst <- split(sr2, idxOdd) # Fwd
    sr1_lst <- sr1_lst[names(sr1_lst) != "NA-NA"]
    
    sr2_lst <- split(sr1, idxOdd) # Rev
    sr2_lst <- sr2_lst[names(sr2_lst) != "NA-NA"]
    
    lapply(seq_along(sr1_lst), function(i) {
      outFile <- file.path(outputDir, names(sr1_lst[i]))
      writeFastq(sr1_lst[[i]], file = paste(outFile, "_R1.fastq.gz", 
                                            sep = ""), mode = mode, compress = T)
      writeFastq(sr2_lst[[i]], file = paste(outFile, "_R2.fastq.gz", 
                                            sep = ""), mode = mode, compress = T)
    })
    
    # Write sequences to FASTQ only undetected reads
    mode <- "w"
    idxNA <- coalesce(gsub("NA-NA", NA, idx), gsub("NA-NA", NA, idxOdd))
    idxNA[is.na(idxNA)] <- "NA-NA"
    
    sumDemultiplex <- c(sumDemultiplex, idxNA[idxNA == "NA-NA"])
    
    sr1_lst <- split(sr2, idxNA) # Fwd
    sr1_lst <- sr1_lst[names(sr1_lst) == "NA-NA"]
    
    sr2_lst <- split(sr1, idxNA) # Rev
    sr2_lst <- sr2_lst[names(sr2_lst) == "NA-NA"]
    lapply(seq_along(sr1_lst), function(i) {
      outFile <- file.path(outputDir, names(sr1_lst[i]))
      writeFastq(sr1_lst[[i]], file = paste(outFile, "_R1.fastq.gz", 
                                            sep = ""), mode = mode, compress = T)
      writeFastq(sr2_lst[[i]], file = paste(outFile, "_R2.fastq.gz", 
                                            sep = ""), mode = mode, compress = T)
    })
    
    progressReport(detail = paste(msg, countReads, "reads done ..."))
  }
  progressReport(detail = paste(msg, "finished, total", countReads, 
                                "demultiplexed reads."))
  close(f1)
  close(f2)
  tab <- as.data.frame(table(sumDemultiplex), stringsAsFactors = F)
  colnames(tab) <- c("BarcodePair", "ReadNumbers")
  outfiles <- list.files(outputDir)
  names(outfiles) <- sub("_R..fastq.gz$", "", outfiles)
  outFileR1 <- outfiles[grep("_R1.fastq.gz$", outfiles)][tab$BarcodePair]
  outFileR2 <- outfiles[grep("_R2.fastq.gz$", outfiles)][tab$BarcodePair]
  tab$FileR1 <- file.path(outputDir, outFileR1)
  tab$FileR2 <- file.path(outputDir, outFileR2)
  return(invisible(tab))
}
