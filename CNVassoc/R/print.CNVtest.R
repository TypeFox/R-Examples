print.CNVtest<-function(x,...){

  if (x$type==1)
    cat("----CNV Wald test----\n")
  else
    cat("----CNV Likelihood Ratio Test----\n")

  cat("Chi=",x$stat,"(df=",x$df,") , pvalue=",x$pvalue,"\n\n")

}