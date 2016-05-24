matplot_enrichment_curve <- function(x, y) {
  if ( dim(x)[1] != length(y) ){
    stop(paste("The number of rows for score matrix must be equal to the number of labels."))
  }
  ncol <- dim(x)[2]
  ### ideal & random ###
  plot( y, y, xlim=c(0, 100), ylim=c(0,100), type="n", las=1,
    xlab="top % of ranked database", ylab="% found Activities (yield)"
  )
  lines( 100 * 0:length(y) / length(y), 100 * c(0, cumsum(sort(y, decreasing=T)) / sum(y)), lwd=2, lty=2)
  lines(0:100, 0:100, lwd=2, lty=3, col="grey")
  
  ### scores ###
  for (j in 1:ncol ){
    ord <- order(x[,j], decreasing=T)
    lines( 100 * 0:length(y) / length(y), 100 * c(0, cumsum(y[ord]) / sum(y)),  lwd=2, col=j+1)
  }
  legend(70, 50, lwd=2, c("ideal", "random",1:ncol), lty=c(2,3,rep(1,length=ncol)), 
    col=c("black", "grey", 2:(ncol+1))
  )
}
