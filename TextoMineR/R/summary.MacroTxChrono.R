summary.MacroTxChrono<-function (object, nword=20, nEig=5, ...) 
{
  res <- object
     if (!inherits(res, "MacroTxChrono")) 
   stop("non convenient object")
    cat("\nTxChrono summary\n")  
    summary(res$Corpus,nword)
    summary(res$res.TxCA,nEig)
    cat("\nCorrelation between chronology and dimensions\n")
    print(res$Correlation)
    ncp=ncol(res$res.TxCA$res.ca$col$coord)
   cat("\n ", ncp, " axes will be taken into account in the constrained hierarchical clustering\n")
    cat("\nHierarchical words\n")
    print(res$HierWord)
  if (!is.null(res$HierSegment)) {
        cat("\nHierarchical segments\n")
        print(res$HierSegment)
 }
    summary(res$VocIndex, nword)    
}
