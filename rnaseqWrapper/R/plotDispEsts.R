plotDispEsts <-
function( cds,cond1,cond2 ){
  plot(
    rowMeans( counts( cds, normalized=TRUE ) ),
    fitInfo(cds)$perGeneDispEsts,
    pch = '.', log="xy", main=paste("Fit of Dispersion Estimate: ",cond1," vs ",cond2,sep=""), 
    xlab="Mean Expression Strength", ylab="Dispersion Values" )
  xg <- 10^seq( -.5, 5, length.out=300 )
  lines( xg, fitInfo(cds)$dispFun( xg ), col="red" )
}
