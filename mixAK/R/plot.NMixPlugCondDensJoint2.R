##
##  PURPOSE:   Plotting of computed pairwise bivariate conditional densities
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   29/05/2009
##  LOG:       08/11/2011:  add.contour, col.add.contour arguments added
##
##  FUNCTIONS: plot.NMixPlugCondDensJoint2 (29/05/2009)
##             
## ======================================================================

## *************************************************************
## plot.NMixPlugCondDensJoint2
## *************************************************************
##
plot.NMixPlugCondDensJoint2 <- function(x, ixcond, imargin, contour=FALSE, add.contour=TRUE, col.add.contour="brown", auto.layout=TRUE, col, lwd=1, main, xylab, ...)
{
  return(plot.NMixPredCondDensJoint2(x=x, ixcond=ixcond, imargin=imargin, contour=contour,
                                     add.contour=add.contour, col.add.contour=col.add.contour,
                                     auto.layout=auto.layout,
                                     col=col, lwd=lwd, main=main, xylab=xylab, ...))

}  
