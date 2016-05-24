##
##  PURPOSE:   Plotting of computed marginal (univariate) densities
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   28/05/2009
##
##  FUNCTIONS: plot.NMixPlugDensMarg (28/05/2009)
##             
## ======================================================================

## *************************************************************
## plot.NMixPlugDensMarg
## *************************************************************
plot.NMixPlugDensMarg <- function(x, auto.layout=TRUE, type="l", col="darkblue", lty=1, lwd=1, main, xlab, ylab, ...)
{
  return(plot.NMixPredDensMarg(x=x, auto.layout=auto.layout, type=type, col=col, lty=lty, lwd=lwd, main=main, xlab=xlab, ylab=ylab, ...))
}

