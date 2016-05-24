plotEcdf <-
function(x, y=NULL, col=c("black","red")) {
  # Graphs one or two empirical cumulative distribution functions.
  # 'x': Vector of numerical observations whose empirical cdf is to be graphed.
  # 'y': Optional vector of observations whose empirical cdf is to be graphed.
  # 'col': Scalar or vector of length 2, specifying the colors of the two empirical distribution functions.
  #    Preferably, the two colors should differ.
  #    \code{col[1]} and \code{col[2]} correspond to \code{x} and \code{y}, respectively.  Type \code{colors()} for selections.
  # example:   plotEcdf( c(2,4,9,6), c(1,7,11,3,8) )
  # example:   plotEcdf( c(2,4,9,6), c(1,7,11,3) )
  if (!is.numeric(x))  stop("'x' must be numeric.")
  if (!is.null(y) & !is.numeric(y))  stop("'y' must be numeric or NULL.")
  if (is.null(y))  plot.ecdf(x, col=col[1])
  if (!is.null(y))  {
     num.grid.points <- 10001;  delta <- 0.1
     xmin <- min(x,y)-(max(x,y)-min(x,y))*delta; xmax <- max(x,y)+(max(x,y)-min(x,y))*delta
     plot.ecdf(x, xlim=c(xmin,xmax), cex=1.3, lwd=5, col=col[1])
     xmin <- min(x,y)-2*(max(x,y)-min(x,y))*delta; xmax <- max(x,y)+2*(max(x,y)-min(x,y))*delta
     f2 <- function(u) { f <- rep(NA, length(u))
             for (i in 1:length(u))  { f[i] <- mean(y<=u[i]) } ;      return(f)  }
     if (length(col)==1)  col <- c(col, col)
     curve(f2, xmin, xmax, num.grid.points, add=TRUE, type="p", pch=20, cex=0.03, col=col[2])
     fy <- function(u){ fy <- mean(y<=u) }
     for (i in 1:length(y))  {
        curve(fy, y[i], y[i], 1, add=TRUE, type="p", pch=20, cex=1.1, col=col[2])  }  }
}
