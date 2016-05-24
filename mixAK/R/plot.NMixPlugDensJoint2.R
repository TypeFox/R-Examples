##
##  PURPOSE:   Plotting of computed pairwise joint densities (plug-in version)
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   28/05/2009
##  LOG:       08/11/2011:  add.contour, col.add.contour arguments added
##
##  FUNCTION:  plot.NMixPlugDensJoint2 (28/05/2009)
##             
## ======================================================================

## *************************************************************
## plot.NMixPlugDensJoint2
## *************************************************************
plot.NMixPlugDensJoint2 <- function(x, contour=FALSE, add.contour=TRUE, col.add.contour="brown", auto.layout=TRUE, col, lwd=1, main, xylab, ...)
{
  return(plot.NMixPredDensJoint2(x=x, K=0, contour=contour,
                                 add.contour=add.contour, col.add.contour=col.add.contour,
                                 auto.layout=auto.layout, col=col, lwd=lwd, main=main, xylab=xylab, ...))
}
