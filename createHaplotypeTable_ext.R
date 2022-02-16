createFinalHaplotypTableExtended <- function(outputDir, sampleTable, markerTable, referenceSequence, snpList, postfix, 
                                     minHaplotypCoverage=3, minReplicate=2, 
                                     detectability=1/100, minSampleCoverage=300,
                                     filterIndel=T){
  # check args
  stopifnot(
    is.character(outputDir), length(outputDir) == 1, file.exists(outputDir),
    is.data.frame(sampleTable), all(c("MarkerID", "ReadFile", "SampleID", "SampleName") %in% colnames(sampleTable)),
    is.character(sampleTable$ReadFile), all(file.exists(sampleTable$ReadFile)),
    is.data.frame(markerTable), all(c("MarkerID") %in% colnames(markerTable)),
    is.null(snpList) | is.list(snpList),
    is.character(postfix), length(postfix) == 1,
    is.numeric(minHaplotypCoverage), length(minHaplotypCoverage) == 1,
    is.numeric(minReplicate), length(minReplicate) == 1,
    is.numeric(detectability), length(detectability) == 1,
    is.numeric(minSampleCoverage), length(minSampleCoverage) == 1
  )
  
  devMode <- getOption("HaplotypR.devel")
  if(is.null(devMode))
    devMode <- T
  outFreqDir <- file.path(outputDir, "frequencyFiles_ext")
  dir.create(outFreqDir, recursive = T)
  res <- lapply(markerTable$MarkerID, function(marker){
    outFreqFiles <- file.path(outFreqDir, marker)
    dir.create(outFreqFiles)
    samTab <- sampleTable[sampleTable$MarkerID == marker,]
    if(marker %in% names(snpList) & !is.null(snpList[[marker]])){
      snpSet <- as.data.frame(snpList[[marker]], stringsAsFactors = FALSE)
      snpSet$Pos <- as.integer(snpSet$Pos)
      snpSet$Alt <- NULL      
    }else{ 
      snpSet <- NULL 
    }
    
    prefix <- sub(".fastq.gz$", "", basename(as.character(samTab$ReadFile)))
    
    
    # Create frequency files and count table
    tab <- createContingencyTable(inputFiles = as.character(samTab$ReadFile), sampleNames=as.character(samTab$SampleID), dereplicated=F,
                                  inputFormat="fastq", outputDir=outFreqFiles, replicatNames=NULL, haplotypeUIDFile=NULL)

    if(devMode)
      write.table(tab, file=file.path(outputDir, sprintf("contingencyTable_%s%s.txt", marker, postfix)), sep="\t", row.names=F, quote=F)
    fnAllSeq <- file.path(outFreqFiles, sprintf("allSequences_%s%s.fasta", marker, postfix))
    file.rename(file.path(outFreqFiles, "allHaplotypes.fa"), fnAllSeq)
    frqfile <- file.path(outFreqFiles, paste(prefix, "_hapFreq.fa", sep=""))
    
    # run cluster with Rswarm package
    outCluster <- file.path(outputDir, "cluster_ext", marker)
    dir.create(outCluster, recursive=T)
    clusterFilenames <- clusterReads(frqfile, outCluster, prefix)
    
    # check for chimeras with Rvsearch package
    chimeraFilenames <- checkChimeras(clusterFilenames[,"RepresentativeFile"], method="vsearch")
    
    # create an overview table
    overviewHap <- createHaplotypOverviewTable(allHaplotypesFilenames=fnAllSeq, clusterFilenames=clusterFilenames, 
                                               chimeraFilenames=chimeraFilenames, referenceSequence=referenceSequence[marker],
                                              snpSet=snpSet, filterIndel=filterIndel)
    if(devMode)
      write.table(overviewHap, file=file.path(outputDir, sprintf("InitialOverviewHap_%s%s.txt", marker, postfix)), sep="\t", row.names=F, quote=F)
    ## create final haplotype
    overviewHap$FinalHaplotype <- NA_character_
    overviewHap[overviewHap$representatives, "FinalHaplotype"] <- marker
    # Set singelton
    overviewHap[overviewHap$singelton, "FinalHaplotype"] <- "Singelton"
    # Set Indel
    overviewHap[!is.na(overviewHap$indels) & overviewHap$indels & overviewHap$representatives, "FinalHaplotype"] <- "Indels"
    # Set chimera
    overviewHap[!is.na(overviewHap$chimeraScore) & is.na(overviewHap$nonchimeraScore), "FinalHaplotype"] <- "Chimera"
    # Cluster identical SNPs pattern
    idx <- overviewHap$FinalHaplotype %in% marker
    snps <- unique(na.omit(overviewHap$snps[idx]))
    if(any(snps=="")){
      hapNames <- paste(marker, 1:sum(idx), sep="-")
      overviewHap$FinalHaplotype[idx] <- paste(marker, 1:sum(idx), sep="-")
    }else if(length(snps) > 0){
      names(snps) <- paste(marker, seq_along(snps), sep="-")
      overviewHap$FinalHaplotype[idx] <- names(snps)[match(overviewHap$snps[idx], snps)]
    }
    overviewHap <- as.data.frame(overviewHap)
    if(devMode)
      write.table(overviewHap, file=file.path(outputDir, sprintf("HaplotypeOverviewTable_%s%s.txt", marker, postfix)), sep="\t")
    
    # getHaplotype sequence
    if(!is.null(referenceSequence[marker]) & filterIndel){
      hapSeq <- lapply(snps, function(sp){
        replaceLetterAt(referenceSequence[[marker]], at=snpSet$Pos, letter=sp)
      })
      hapSeq <- DNAStringSet(hapSeq)
      writeFasta(hapSeq, file.path(outputDir, file=sprintf("%s_HaplotypeSeq%s.fasta", marker, postfix)))
    } else {
      idx <- grep(marker, overviewHap$FinalHaplotype)
      hapSeq <- readFasta(file.path(outputDir, "frequencyFiles", marker, file=sprintf("allSequences_%s%s.fasta", marker, postfix)))
      hapSeq <- hapSeq[vctrs::vec_group_id(hapSeq) %in% overviewHap$HaplotypesName[idx]]
      hapSeq <- sread(hapSeq)
      names(hapSeq)  <- overviewHap$FinalHaplotype[idx]
      writeFasta(hapSeq, file.path(outputDir, file=sprintf("%s_HaplotypeSeq%s.fasta", marker, postfix)))
    }
    repfile <- clusterFilenames[,"RepresentativeFile"]
    
    hab <- rep(0, dim(overviewHap)[1])
    names(hab) <- rownames(overviewHap)
    
    haplotypesSample <- lapply(seq_along(repfile), function(i){
      sr1 <- readFasta(repfile[i])
      vec <- do.call(rbind, strsplit(as.character(ShortRead::id(sr1)), "_"))
      clusterResp <- vec[,1]
      clusterSize <- as.integer(vec[,2])
      hab[clusterResp] <- clusterSize
      hab
    })
    haplotypesSample <- do.call(cbind, haplotypesSample)
    haplotypesSample <- haplotypesSample[rowSums(haplotypesSample)>0, , drop = FALSE]
    colnames(haplotypesSample) <- sub(".representatives.fasta", "", basename(repfile))
    
    overviewHap <- overviewHap[rownames(haplotypesSample),]
    
    # set final haplotype names
    rownames(haplotypesSample) <- overviewHap[rownames(haplotypesSample), "FinalHaplotype"]
    if(devMode)
      write.table(cbind(HaplotypNames=rownames(haplotypesSample),haplotypesSample), 
                  file=file.path(outputDir, sprintf("rawHaplotypeTable_%s%s.txt", marker, postfix)), sep="\t", row.names=F, col.names=T)
    
    haplotypesSample <- reclusterHaplotypesTable(haplotypesSample)
    if(devMode)
      write.table(cbind(HaplotypNames=rownames(haplotypesSample),haplotypesSample), 
                  file=file.path(outputDir, sprintf("reclusteredHaplotypeTable_%s%s.txt", marker, postfix)), sep="\t", row.names=F, col.names=T)
    
    # Apply cut-off haplotype only in 1 sample
    haplotypesSample <- callHaplotype(haplotypesSample, minHaplotypCoverage=minHaplotypCoverage, minReplicate=minReplicate, 
                                      detectability=detectability, minSampleCoverage=1)
    
    if(devMode) 
      write.table(cbind(HaplotypNames=rownames(haplotypesSample), haplotypesSample), 
                  file=file.path(outputDir, sprintf("finalHaplotypeTable_Hcov%.0f_Scov%.0f_occ%i_sens%.4f_%s%s.txt",
                                                    minHaplotypCoverage, minSampleCoverage, minReplicate, detectability, marker, postfix)),
                  sep="\t", row.names=F, col.names=T)
    
    
    # check replicates
    idx <- split(1:dim(haplotypesSample)[2], samTab$SampleID)
    markerRes <- lapply(idx, function(i){
      tab <- callHaplotype(haplotypesSample[,i, drop=F], minHaplotypCoverage=minHaplotypCoverage, 
                           minReplicate=minReplicate, detectability=detectability, minSampleCoverage=minSampleCoverage, 
                           reportBackground=T)
      tab <- cbind(samTab[rep(i,each=dim(tab)[1]), c("SampleID","SampleName","MarkerID")], 
                   Haplotype=rownames(tab), Reads=as.integer(tab), FlagChimera=F)
      colnames(tab) <- c("SampleID","SampleName","MarkerID","Haplotype","Reads")
      rownames(tab) <- NULL
      
      #check individual samples for chimera
      do.call(rbind, lapply(split(tab, tab$SampleID), function(tt){
        chim <- NULL
        hIdx <- grep(marker, tt$Haplotype)
        if(length(hIdx)>2){
          chim <- flagChimera(tt[hIdx,], overviewHap)
        }
        tt$FlagChimera <- tt$Haplotype %in% chim
        return(tt)
      }))
      return(tab)
    })
    markerRes <- do.call(rbind.data.frame, markerRes)
    rownames(markerRes) <- NULL
    markerResFn <- file.path(outputDir, sprintf("finalHaplotypeList_Hcov%.0f_Scov%.0f_occ%i_sens%.4f_%s%s.txt", 
                                                minHaplotypCoverage, minSampleCoverage, minReplicate, detectability, marker, postfix))
    write.table(markerRes, file=markerResFn, sep="\t", row.names=F, col.names=T)
    return(markerRes)
  })
  names(res) <- markerTable$MarkerID
  return(res)
}
