################################################################################
#  Plot of HRs' with confidence intervals                                      #
################################################################################
#                                                                              #
#  Plots the hazard ratios' estimated values and confidence intervals for      #
#  parametric frailty models                                                   #
#                                                                              #
#  Its parameters are                                                          #
#   - x        : the frailty model, object of class 'predict.parfm'            #
#   - xlim     : xlim (see pars())                                             #
#   - main     : main (see pars())                                             #
#   - coef     : the name of the covariates whose coefficients must be plotted #
#   - signcol  : Should non-significant coefficients be plotted in grey?       #
#   - cex      : cex (see pars())                                              #
#                                                                              #
#   Date: February 21, 2012                                                    #
#   Last modification on: October 17, 2012                                     #
################################################################################

plot.parfm <- function(x, 
                       xlim=NULL,
                       main=NULL,
                       coef=NULL,
                       signcol=TRUE,
                       cex=1,
                       ...){
  if (is.null(coef)) {
    n <- sum(!is.na(x[, "p-val"]))
    coef <- rownames(x)[nrow(x) + 1 - n:1]
  } else
    n <- length(coef)
  coef <- coef[n:1]
  
  intervals <- cbind(exp(x[coef, "ESTIMATE"]),
                     ci.parfm(x)[coef, , drop=FALSE])
  
  if (is.null(xlim)) {
    range <- c(0, 
               max(max(intervals[is.finite(intervals)], na.rm=TRUE) * 1.2, 2))
    if (!is.finite(range[2])) 
      range <- c(0, 100)
  } else {
    range <- xlim
  }
  
  if (is.null(main)) {
    frailty <- attr(x, "frailty")
    dist <- attr(x, "dist")
    main <- paste(c(none="Cox ", 
                    paste(c("Gamma", "Inverse Gaussian", "Positive stable",
                            "Lognormal"),
                  "frailty "))[which(c("none", "gamma", "ingau", "possta",
                                       "lognormal") == frailty)], 
                  "model\nwith ", 
                  dist, " baseline", sep="", collapse="")
  }
  
  if (signcol) {
    color <- apply(intervals[, c("low", "up"), drop=FALSE], 1, 
                   function(int) {
                     if (is.na(int[1]) | is.na(int[2]) | (int[1]<1 & int[2]>1))  
                       "grey" else "black"
                   })
  } else 
    color <- "black"
  
  par(mar=c(5, 5, 4, 2) + .1)
  
  plot(as.vector(intervals), rep(1:n, 3),
       xlim=range,
       ylim=c(.2, n + .8),
       bty="]",
       xlab="HR", 
       ylab="", 
       yaxt="n",
       main=main,
       pch= c(rep(20, n), rep(91, n), rep(93, n)),
       cex=cex,
       col=color)
  
  abline(h=1:n, v=0:1, 
         lty=c(rep(3, n), 1, 2), 
         col=c(rep("grey", n), rep("black", 2)))
  
  intervals[!is.finite(intervals[,3]),3] <- 
    intervals[!is.finite(intervals[,3]),1] + 10^10
  
  segments(intervals[, 2], 1:n, intervals[, 3], 1:n, col=color)
  
  mtext(substr(dimnames(intervals)[[1]], 1, 10),
        side=2,
        at=1:n,
        las=1,
        cex=.7,
        col=color)	
}
