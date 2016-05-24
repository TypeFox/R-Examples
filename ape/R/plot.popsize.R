## plot.popsize.R (2004-07-4)

##   Plot population size in dependence of time

## Copyright 2004 Rainer Opgen-Rhein and Korbinian Strimmer

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

plot.popsize <- function(x, show.median=TRUE,
    show.years=FALSE, subst.rate, present.year, ...)
{
  if (class(x) != "popsize")
    stop("object \"x\" is not of class \"popsize\"")

  ylim <- c(min(popsize[,2:5]),max(popsize[,2:5]))
  if (show.years)
  {
    x1 <- -x[,1]/subst.rate+present.year
    xlab <- "time (years)"
    xlim <- c(min(x1),max(x1))
  }
  else
  {
    x1 <- x[,1]
    xlab <- "time (past to present in units of substitutions)"
    xlim <- c(max(x1),min(x1))
  }

  if (show.median)
    plot(x1,x[,3],type="s", xlim=xlim, ylim=ylim, xlab=xlab,ylab="effective population size",log="y", lwd=2.5, ...) #median
  else
    plot(x1,x[,2],type="s", xlim=xlim, ylim=ylim, xlab=xlab,ylab="effective population size",log="y", lwd=2.5, ...) #median

  lines(x1,x[,4], ...)
  lines(x1,x[,5], ...)
}



lines.popsize <- function(x, show.median=TRUE,
    show.years=FALSE, subst.rate, present.year, ...)
{
  if (class(x) != "popsize")
    stop("object \"x\" is not of class \"popsize\"")

  if (show.years)
  {
    x1 <- -x[,1]/subst.rate+present.year
  }
  else
  {
    x1 <- x[,1]
  }


  if (show.median)
    lines(x1,x[,3], lwd=2.5, ...) #median
  else
    lines(x1,x[,2], lwd=2.5, ...) #median

  lines(x1,x[,4], ...)
  lines(x1,x[,5], ...)
}
