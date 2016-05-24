rootonorm <- function(x, breaks="Sturges",
                      type=c("hanging", "deviation"),
                      scale=c("sqrt", "raw"),
                      zeroline=TRUE,                      
                      linecol="red", rectcol="lightgrey",
                      xlab=xname,
                      ylab="Sqrt(frequency)",
                      yaxt="n",                      
                      ylim=NULL,
                      mu=mean(x), s=sd(x),
                      gap=0.1, ...) {


  if (!is.numeric(x)) 
    stop("'x' must be numeric")

  # Fix the xlabel if it isn't specified
  xname <- deparse(substitute(x))
 
  scale <- match.arg(scale)
  if (is.character(scale) && scale == "raw") {
    scale <- match.fun("as.numeric")
    if (missing(ylab)) {
      ylab <- "Frequency"
    }
  }
  else {
    scale <- match.fun(scale)
  } 

  type <- match.arg(type)
  
  h <- hist(x, breaks=breaks, plot=FALSE)
  if (!h$equidist) stop("breaks must be equally spaced")
  
  nbins <- length(h$counts)
  nobs <- sum(h$counts)

  expected <- nobs*diff(pnorm(h$breaks, mu, s))
  
  d.gap <- min(diff(h$breaks)) * gap /2
 
  plot.range <- range(pretty(h$breaks))  
  z <- seq(plot.range[1], plot.range[2], length.out=200)
  z.y <- min(diff(h$breaks))*nobs*dnorm(z, mu, s)
  
  minval <- min(scale(expected)-scale(h$counts))

  if (is.null(ylim)) {
    ylim <- c(minval, scale(max(expected,z.y)))
  }
  
  plot(z, z, type="n",
       xlab=xlab,
       ylab=ylab,
       yaxt=yaxt,
       ylim=ylim,
       ...)

  if (type=="deviation") {  
    for(i in 1:nbins) {
      rect(h$breaks[i]+d.gap, scale(expected[i])-scale(h$counts[i]),
           h$breaks[i+1]-d.gap, 0, col=rectcol)
    }    
  }
  else {
    for(i in 1:nbins) {
      rect(h$breaks[i]+d.gap, scale(expected[i])-scale(h$counts[i]),
           h$breaks[i+1]-d.gap, scale(expected[i]), col=rectcol)
    }
  }
  
  lines(z, scale(z.y), col=linecol, ...)
  if (zeroline) {
    abline(h=0, lty=3)
  }

  invisible(h$counts)
}
