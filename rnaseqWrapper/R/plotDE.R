plotDE <-
  function( res,cond1,cond2 ){
  plot(res$baseMean,
       res$log2FoldChange,
       log="x", pch=20, cex=.3,
       main= paste("Log2 Fold Change Versus Base Means (5% FDR): ",cond1," vs ",cond2,sep=""), sub="MA-Plot",
       xlab="Base Means", ylab="Log2 Fold Change",
       col = ifelse( res$padj < .05, "red", "black" ) )
}
