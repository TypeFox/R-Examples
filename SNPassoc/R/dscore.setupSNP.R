dscore.setupSNP <- function(x, ...)
 {
   ss <- summary(x, print=FALSE)
   MAF <- 1 - (ss$major.allele.freq)/100
   ans <- dscore(MAF)
   attr(ans, "probs") <- MAF
   ans   
 }
