##
##  PURPOSE:   Plotting of computed univariate conditional densities
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   28/05/2009
##
##  FUNCTIONS: plot.NMixPlugCondDensMarg (28/05/2009)
##             
## ======================================================================

## *************************************************************
## plot.NMixPlugCondDensMarg
## *************************************************************
##
plot.NMixPlugCondDensMarg <- function(x, ixcond, imargin, over=FALSE, auto.layout=TRUE, type="l", lwd=1, lty, col, main, xlab, ylab, ylim, annot=TRUE, ...)
{
  return(plot.NMixPredCondDensMarg(x=x, ixcond=ixcond, imargin=imargin, over=over, auto.layout=auto.layout, type=type,
                                   lwd=lwd, lty=lty, col=col, main=main, xlab=xlab, ylab=ylab, ylim=ylim, annot=annot, ...))  
}

