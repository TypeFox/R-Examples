`DistributionPlotBinomial` <-
function(size=200,prob=0.5, xlab="Number of Successes",
           ylab="Probability Mass", signif.digits=3,
           main=paste("Binomial Distribution: n =",size,"p =",
             signif(prob,digits=signif.digits)) )
{
  ## Purpose: plot binomial distribution (guts borrowed from Rcmdr)
  ##
  ## size:     number of trials
  ## prob:     probability of success
  ## main:     title for plot
  ## ylab:     y-axis label
  ## xlab:     x-axis label
  ## signif.digits:   number of significant digits for default 'main' title

  min <- qbinom(.0005, size=size, prob=prob)
  max <- qbinom(.9995, size=size, prob=prob)
  .x <- min:max
  plot(.x, dbinom(.x, size=size, prob=prob), xlab=xlab, ylab=ylab,
       main=main, type="h")
  points(.x, dbinom(.x, size=size, prob=prob), pch=16)
  abline(h=0, col="gray")
}

