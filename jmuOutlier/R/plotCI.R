plotCI <-
function(CI, mu=NULL, plot.midpoints=TRUE,
                    col=c("black", "red", "darkgreen", "purple"))  {
  # Plots a set of confidence intervals on the same graph.
  # 'CI':  N by 2 matrix or 2 by N matrix consisting of N two-sided confidence intervals.
  # 'mu':  Value of the estimated parameter, and is the population mean for t-confidence intervals.
  # 'plot.midtpoints': Logical; plots the midpoints of the confidence intervals if \code{TRUE} (default);
  #                             otherwise, does not plot the midpoints.
  # 'col': Scalar or vector, specifying the colors of the
  #    plotted confidence intervals, the vertical line going through \code{mu}, and the midpoints.
  #    Type \code{colors()} for selections.
  # Note:  Confidence intervals not containing \code{mu} are in the color \code{col[2]}.
  #        Confidence intervals containing \code{mu} are in the color \code{col[3]}.
  # example:  # Plot 20 confidence intervals, each based on 13 observations from a Normal(mean=70,sd=10) distribution.
  #           CI = replicate( 20, CI.t.test( rnorm( 13, 70, 10 ) ) )
  #           plotCI( CI, 70 )
  if (!is.numeric(CI))  stop("'CI' must be numeric.")
  if (length(col)==1) {col <- rep(col, 4)};  if (length(col)==2) {col <- c(col, col)}
  if (length(col)==3) {col <- c(col, col[1])}
  if (is.vector(CI)) { CI <- t(matrix(CI)) }
  if (dim(CI)[1]==2 & dim(CI)[2]!=2)   CI <- t(CI)
  if (dim(CI)[2]!=2)  stop("'CI' must be an N by 2 matrix of N two-sided confidence intervals.")
  if ( (!is.null(mu) & !is.numeric(mu)) | length(mu)>1 )  stop("'mu' must be a scalar numeric or NULL.")
  if (!is.logical(plot.midpoints))  stop("'plot.midpoints' must be logical.")
  main <- ifelse( is.null(mu), "", paste(c( prettyNum(100*mean( (CI[,1]<=mu)*(mu <=CI[,2]) )),
                   "% of these confidence intervals contain ",prettyNum(mu),"." ), collapse =" ")   )
  plot(NA, NA, xlim=c(min(CI[,1],mu),max(CI[,2],mu)), ylim=c(0.5,dim(CI)[1]+0.5) ,
       type="s", main=main,
       xlab=paste(c("CONFIDENCE INTERVAL",ifelse(dim(CI)[1]==1,"","S")),collapse =""),
       ylab="label")
  if (is.null(mu))  {col1 <- rep(col[3], dim(CI)[1])}
  else {col1 <- ifelse(CI[,1]<=mu & mu<=CI[,2], col[3], col[2])}
  for (i in 1:dim(CI)[1]) {lines( c(CI[i,1], CI[i,2]), c(i,i), col=col1[i])}
  midpoint <- (CI[,1]+CI[,2])/2 ;   delta <- (max(CI[,2],mu)-min(CI[,1],mu))*0.002
  if (!is.null(mu))  { lines( c(mu,mu), c(0.5,dim(CI)[1]+0.5), col=col[1])   }
  if (plot.midpoints) { for (i in 1:dim(CI)[1]) {lines( c(midpoint[i]-delta,midpoint[i]+delta), 
                                                        c(i,i), lwd=4, col=col[4])} }
}
