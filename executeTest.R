library("HaplotypR")
library("ShortRead")

# Define output directory 
outputDir <- "test"  
# Create output directory
if(!dir.exists(outputDir))
  dir.create(outputDir, recursive=T)

# Copy example files to output directory
file.copy(from=system.file(package="HaplotypR", "extdata"), to="./test", recursive = T)

# List files example files in output direcoty
dir(file.path("extdata"))


############## Run demultiplexing by sample and rename output files
# set input file path
primerFile <- "/home/rstudio/markerFile.txt"
sampleFile <- "/home/rstudio/test/test_sampleFile.txt"
fnBarcodeF <- "/home/rstudio/barcode_Fwd.fasta"
fnBarcodeR <- "/home/rstudio/barcode_Rev.fasta"
reads <- list.files("/home/rstudio/test/", pattern="Barcodes", full.names = T)

# source new funs
source("./haplotypR_ext/deplexBySample_ext.R")
source("./haplotypR_ext/extendSearch_ext.R")


# create output subdirectory 
outDeplexSample <- file.path(outputDir, "dePlexSampleTest")
dir.create(outDeplexSample)

# None
dePlexSample <- deplexSampleExtended(reads[1], reads[2], fnBarcodeF, fnBarcodeR, outDeplexSample, max.mismatch = 0, with.indels = F)

# All
outDeplexSample <- file.path(outputDir, "dePlexSampleTestIndels")
dir.create(outDeplexSample)
dePlexSamplewIndels <- deplexSampleExtended(reads[1], reads[2], fnBarcodeF, fnBarcodeR, outDeplexSample, max.mismatch = 1, with.indels = T)

# All SNPs but only deletion 1st, 7th, 8th
outDeplexSample <- file.path(outputDir, "dePlexSampleTestMM1")
dir.create(outDeplexSample)
dePlexSamplewOneMM <- deplexSampleExtended(reads[1], reads[2], fnBarcodeF, fnBarcodeR, outDeplexSample, max.mismatch = 1, with.indels = F)
