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
 
rapidr.dir <- system.file(package = "RAPIDR")
gcContent.fname <- paste(rapidr.dir, "/data/gcContent.RData", sep = "") 
ref.set  <- createReferenceSetFromCounts(counts.fname, outcomes, gcCorrect = FALSE, PCA = FALSE, method = "NCV") 
#ref.set  <- createReferenceSetFromCounts(counts.fname, outcomes, gcCorrect = TRUE, PCA = FALSE, method = "zscore", gcContentFile = gcContent.fname)
#unknowns.results <- testUnknowns(ref.set, counts.fname, gcContentFile = gcContent.fname)
unknowns.results <- testUnknowns(ref.set, counts.fname) 
print(unknowns.results[['results']])  
 

