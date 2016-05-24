plot.meanssaemodel <-
function(x, y=NULL, type="e", sort=NULL, ...){
   # y is part of the generic (later, we will use y to allow plot comparison)
   # extract robustness tuning constant
   k <- attr(x, "robustness")
   fe <- x$fixeff
   re <- x$raneff
   means <- x$means
   g <- length(re)
   areaNames <- rownames(re)
   mspe <- x$mspe
   # sorting
   if (!is.null(sort)){
      ord <- switch(sort, raneff=order(re), ranef=order(re), fixeff=order(fe), fixef=order(fe), means=order(means)) 
      re <- re[ord]
      fe <- fe[ord]
      means <- means[ord]
      areaNames <- areaNames[ord]
      if (!is.null(mspe)){
	 mspe <- mspe[ord]
      }
   }
   # prepare the plot
   at <- 1:g
   # without mspe
   if (is.null(mspe)){
      ra <- range(means)
      ra[1] <- min(fe) 
      # add an null-line to the plot (for the legend) 
      g <- g + 1
      re <- c(re, NA)
      fe <- c(fe, NA)
      means <- c(means, NA)  
      areaNames <- c(areaNames, "") 
      op <- par(mfcol=c(1,1), mar=c(4,8,2,4))
      plot(means, at, type="b", col=2, lwd=2, axes=FALSE, xlab="predicted mean", ylab="", xlim=ra, main="Predicted means")
      lines(fe, at, type="b", col=1, lwd=2, xlim=ra)
      axis(2, seq(1, g), labels=areaNames, las=1)
      axis(1)
      grid(col="gray65", lty=2, lwd=1)
      box()
      legend("top", pch=c(1,1), lty=c(1,1), col=c(1,2), legend=c("fixeff prediction", "full prediction"), bg="white", ncol=2) 
   }else{
      # with mspe
      # check type
      if(is.na(match(type, c("e", "l")))) stop("plot type must be either 'e' or 'l' \n")
      # compute plotting range
      ra <- range(means+sqrt(mspe), means-sqrt(mspe))
      op <- par(mfcol=c(1,1), mar=c(8,4,2,4))
      plot(at, means, type="b", col=2, lwd=2, axes=FALSE, xlab="", ylab="predicted area mean", ylim=ra, main="Predicted means (+/- SQRT[MSPE])")
      if (type == "l"){
	 lines(at, means+sqrt(mspe), type="b", col=1, lwd=2, lty=2, xlim=ra)
	 lines(at, means-sqrt(mspe), type="b", col=1, lwd=2, lty=2, xlim=ra)
      }
      if (type == "e"){
	 arrows(at, means-sqrt(mspe), at, means+sqrt(mspe), angle=90, code=3, ...)
      }
      axis(1, seq(1, g), labels=abbreviate(areaNames, minlength=12), las=2)
      axis(2)
      grid(col="gray65", lty=2, lwd=1)
      box()
   }
}

