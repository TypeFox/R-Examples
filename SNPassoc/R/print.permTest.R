print.permTest<-function(x, level=0.05, digits=8, ...)
 {
  if (!inherits(x, "permTest")) 
      stop("x must be an object of class 'permTest'")

  cat("\n")
  cat(paste("Permutation test analysis ","(",(1-level)*100,"% confidence level)",sep=""), "\n")
  cat("------------------------------------------------ \n")
  cat("Number of SNPs analyzed: ",sum(x$nSNPs),"\n")
  cat("Number of valid SNPs (e.g., non-Monomorphic and passing calling rate): ",x$nSNPs[1],"\n")
  cat("P value after Bonferroni correction: ",round(level/x$nSNPs[1] ,digits),"\n")

  control<-attr(x,"method")
  if (control==1)
   {
    pos<-ceiling(length(x$pmin)*(1-level))
    pPerm<-sort(x$pmin,decreasing=TRUE)[pos]
    cat("\n")
    cat("P values based on permutation procedure: \n")
    cat("P value from empirical distribution of minimum p values: ",round(pPerm ,digits),"\n")
    cat("P value assuming a Beta distribution for minimum p values: ",round(x$psig ,digits),"\n")
   }
  if (control==2)
   {
    cat("\n")
    cat(paste("Rank truncated product of the K=", x$K, " most significant p-values:",sep=""))
    cat("\n")
    cat("Product of K p-values (-log scale): ",x$rtp ,"\n")
    cat("Significance: ",ifelse(x$sig==0,"<0.001",x$sig) ,"\n")
   }
 
 }