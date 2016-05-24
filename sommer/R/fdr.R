fdr <- function(p, fdr.level=0.05){
  ##### transform to real p-values
  pval <- 10^-p
  ##### adjust for FDR
  pvalA <- p.adjust(pval, method="fdr")
  #plot(pvalA)
  ##### transformed brough back to original scale
  pvalA.l10 <- -log(pvalA, 10)
  #plot(pvalA.l10)
  ##### for the transformed p-values (-log10) we can use an easy p-val
  fdr.ad <- -log(fdr.level, 10)
  ## FDR in original scale
  vv <- which(pvalA.l10 > fdr.ad)
  if(length(vv)>0){
    fdr.or <- min(p[vv])
  }else{
    fdr.or <- NULL
  }
  ######
  result <- list(p.ad=pvalA.l10, fdr=fdr.ad, p.or=p, fdr.or=fdr.or )
  return(result)
}