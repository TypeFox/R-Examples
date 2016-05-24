recode <- function (X, alleles, values = c(0, 1, 2), pos1 = 1, pos2 = 3, minor = FALSE) 
{
  if(is.matrix(X) | is.data.frame(X)) {
    n <- nrow(X)
    p <- ncol(X)
    Y <- NULL
    for (i in 1:p) {
      snp <- rep(NA, n)
      al1 <- substr(alleles[i], pos1, pos1)
      al2 <- substr(alleles[i], pos2, pos2)
      if (i%%100 == 0) 
        cat("Converting marker ", i, "\n")
      hom1 <- paste(al1, al1, sep = "")
      hom2 <- paste(al2, al2, sep = "")
      het1 <- paste(al1, al2, sep = "")
      het2 <- paste(al2, al1, sep = "")
    
      snp[X[, i] == hom1] <- values[1]
      snp[X[, i] == het1 | X[, i] == het2] <- values[2]
      snp[X[, i] == hom2] <- values[3]
    
      if(minor) {# coding 0,1,2 reflects the number of minor alleles.
         nhom1 <- sum(X[, i] == hom1) # count the numbers of each homozygote
         nhom2 <- sum(X[, i] == hom2)
         al1minor <- nhom1 <= nhom2 # first allele is minor allele?
         if(al1minor) {
           snp[X[, i] == hom1] <- values[3]
           snp[X[, i] == hom2] <- values[1]
         } else {
           snp[X[, i] == hom1] <- values[1]
           snp[X[, i] == hom2] <- values[3]
         }
      }
    Y <- cbind(Y, snp)
  }
  colnames(Y) <- paste("SNP", 1:p, sep = "")
  return(Y)
  } # if(is.matrix)...
  if(is.vector(X)) {
     n <- length(X)
     snp <- rep(NA, n)
     al1 <- substr(alleles, pos1, pos1)
     al2 <- substr(alleles, pos2, pos2)
     hom1 <- paste(al1, al1, sep = "")
     hom2 <- paste(al2, al2, sep = "")
     het1 <- paste(al1, al2, sep = "")
     het2 <- paste(al2, al1, sep = "")
  
     snp[X == hom1] <- values[1]
     snp[X == het1 | X == het2] <- values[2]
     snp[X == hom2] <- values[3]
  
     if(minor) {# coding 0,1,2 reflects the number of minor alleles.
       nhom1 <- sum(X == hom1) # count the numbers of each homozygote
       nhom2 <- sum(X == hom2)
       al1minor <- nhom1 <= nhom2 # first allele is minor allele?
       if(al1minor) {
         snp[X == hom1] <- values[3]
         snp[X == hom2] <- values[1]
       } else {
         snp[X == hom1] <- values[1]
         snp[X == hom2] <- values[3]
       }
     } # end if minor
     return(snp)
  } # end if vector
  if(!is.matrix(X) & !is.data.frame(X) & !is.vector(X)) {
    stop("recode: input is not a matrix, a dataframe or a vector\n")
  }
}
