plot.SimCi <-
function (x, xlim, xlab, ylim, ...) {


comparison <- rep(x$comp.names, each=length(x$resp))
endpoint <- rep(x$resp, times=length(x$comp.names))
estimate <- lower.raw <- upper.raw <- lower <- upper <- NULL
for (i in 1:length(x$comp.names)) {
  estimate <- c(estimate, x$estimate[i, ])   
  lower <- c(lower, x$lower[i, ])
  upper <- c(upper, x$upper[i, ])
}

out <- data.frame(comparison, endpoint, estimate, lower, upper)
opar <- par(no.readonly=TRUE)                                       # make a copy of current settings    
for (i in 1:length(x$resp)) {   
  if (i > 1) par(ask=TRUE)   
  r <- x$resp[i]
  outEP <- out[out$endpoint==r,]

  xrange <- c(min(outEP$lower), max(outEP$upper))       
  if (!is.finite(xrange[1])) xrange[1] <- min(outEP$estimate)
  if (!is.finite(xrange[2])) xrange[2] <- max(outEP$estimate)   
  yvals <- length(outEP$comparison):1
  xlim <- xrange    
  ylim <- c(0.5, length(outEP$comparison) + 0.5)           
  plot(c(outEP$lower, outEP$upper), rep.int(yvals, 2), type="n", 
       axes=FALSE, xlab="", ylab="", xlim=xlim, ylim=ylim, ...)
  axis(1, ...)
  axis(2, at=length(outEP$comparison):1, labels=outEP$comparison, las=1, ...)
  abline(h=yvals, lty=1, lwd=1, col="lightgray")
  vl <- if ("NSD" %in% attributes(x)$names) {1} else {0}          # SimCiRat has attr "NSD"
  abline(v=vl, lty=2, lwd=1, ...)    
  left  <- outEP$lower; left[!is.finite(left)]   <- min(c(0, xrange[1] * 2))
  right <- outEP$upper; right[!is.finite(right)] <- max(c(0, xrange[2] * 2))
  segments(left, yvals, right, yvals, ...)        
  points(outEP$lower, yvals, pch="(", ...)
  points(outEP$upper, yvals, pch=")", ...)
  points(outEP$estimate, yvals, pch=20, ...)      
  main <- list(...)$main
  if (is.null(main)) main <- paste0(format(100 * x$conf.level, 2),
                                   "% confidence level, ",r)     
  if (missing(xlab)) 
    xlab <- if ("NSD" %in% attributes(x)$names) {
      "Ratio"
    } else { 
      "Linear Function"
    } 
  title(main=main, xlab=xlab)
  box()
}

par(opar)                                                           # restore original settings


}
