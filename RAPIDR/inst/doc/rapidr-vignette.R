### R code from vignette source 'rapidr-vignette.rnw'

###################################################
### code chunk number 1: makeBinnedCountsFile (eval = FALSE)
###################################################
## makeBinnedCountsFile(bam.file.list = c("file1.bam", "file2.bam"),
##                      sampleIDs = c("sample1", "sample2"),
##    	             binned.counts.fname = output.fname, 
##                      k = 20000) 


###################################################
### code chunk number 2: makeBinnedCountsFile (eval = FALSE)
###################################################
## mask <- "mask.bed"
## makeBinnedCountsFile(bam.file.list = c("file1.bam", "file2.bam"),
##                      sampleIDs = c("sample1", "sample2"),
##    	             binned.counts.fname = output.fname, 
##                      k = 20000, mask) 


###################################################
### code chunk number 3: createReferenceSetFromCounts
###################################################
library(RAPIDR)
rapidr.dir <- system.file(package = "RAPIDR")
data(outcomes)
data(gcContent) 
# Make some binned data 
T21.pos <- which(outcomes$Dx == "T21")
chr.lens <- sapply(gcContent, length)
chr.names <- names(chr.lens)

# Make the header 
header <- c("SampleID")
for (i in 1:length(chr.lens)) {
   header <- c(header, rep(chr.names[i], chr.lens[i]))
}
nbins <- sum(chr.lens)
ncols <- nbins + 1

binned.counts <- matrix(nrow = nrow(outcomes), ncol = ncols)
for (i in 1:nrow(binned.counts)) {
   binned.counts[i,] <- rpois(ncols, lambda = 100)
   if (i %in% T21.pos) {
      binned.counts[i, 139087:141493] <- rpois(chr.lens[21], lambda = 115)
   }
}
binned.counts[,1] <- outcomes$SampleID
colnames(binned.counts) <- header
t <- tempfile()
write.table(binned.counts, file = t, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ",")
counts.fname <- t
message(t)  
gcContent.fname <- paste(rapidr.dir, "/data/gcContent.RData", sep = "")
head(outcomes)
ref.set <- createReferenceSetFromCounts(counts.fname,
                             outcomes,
   	                         gcCorrect = FALSE,
   	                         PCA = FALSE,
                                 filterBin = FALSE, 
   	                         gcContentFile = gcContent.fname) 


###################################################
### code chunk number 4: print
###################################################
print(ref.set)  


###################################################
### code chunk number 5: testUnknowns
###################################################
test.results <- testUnknowns(ref.set,
                             counts.fname,
                             gcContentFile = gcContent.fname)  
print(test.results)


###################################################
### code chunk number 6: plotTestSample
###################################################
plotTestSample(test.results, "1002")


###################################################
### code chunk number 7: session
###################################################
sessionInfo()


