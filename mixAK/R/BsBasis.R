##
##  PURPOSE:   Create a B-spline basis for specific dataset
##             (with common boundary knots and equidistant inner knots)
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   08/07/2008 (as a stand alone function)
##             26/11/2010:  added to the mixAK package
##
##  FUNCTIONS: BsBasis
##
## ==========================================================================

## *************************************************************
## BsBasis
## *************************************************************
##
BsBasis <- function(degree, ninner, knotsBound, knots, intercept=FALSE, x, tgrid,
                    Bname="B", plot=FALSE, lwd=1, col="blue", xlab="Time", ylab="B(t)",
                    pch=16, cex.pch=1, knotcol="red")
{
  ## Knots
  if (missing(knots)){
    dist <- (knotsBound[2] - knotsBound[1])/(ninner+1)
    knots <- seq(knotsBound[1], knotsBound[2], length=ninner+2)
    knotsInner <- knots[-c(1, length(knots))]
  }else{
    knotsBound <- knots[c(1, length(knots))]
    knotsInner <- knots[-c(1, length(knots))]
  }  

  ## B-spline basis
  BsData <- splines::bs(x, degree=degree, knots=knotsInner, Boundary.knots=knotsBound, intercept=intercept)
  nBs <- ncol(BsData)
  labelBs <- paste(Bname, 1:nBs, sep="")

  attr(BsData, "knots") <- NULL
  attr(BsData, "Boundary.knots") <- NULL
  class(BsData) <- "matrix"
  colnames(BsData) <- labelBs

  ## Figure of the B-spline basis
  if (missing(tgrid)) tgrid <- seq(knotsBound[1], knotsBound[2], length=1000)
  Bs <- splines::bs(x=tgrid, degree=degree, knots=knotsInner, Boundary.knots=knotsBound, intercept=intercept)
  
  if (plot){
    #par(mfrow=c(1, 1), bty="n")
    matplot(x=tgrid, y=Bs, type="l", col=col, lwd=lwd, xaxt="n", xlab=xlab, ylab=ylab)
    axis(1, at=round(c(knotsBound[1], attr(Bs, "knots"), knotsBound[2]), 3))
    points(attr(Bs, "knots"), rep(0, length(attr(Bs, "knots"))), col=knotcol, pch=pch, cex=cex.pch)
    points(knotsBound, rep(0, 2), col=knotcol, pch=pch, cex=cex.pch)
  }  

  attr(Bs, "degree") <- NULL
  attr(Bs, "knots") <- NULL
  attr(Bs, "Boundary.knots") <- NULL
  attr(Bs, "intercept") <- NULL  
  class(Bs) <- "matrix"
  colnames(Bs) <- labelBs
  
  attr(BsData, "knots") <- knots   
  attr(BsData, "knotsInner") <- knotsInner  
  attr(BsData, "knotsBound") <- knotsBound
  attr(BsData, "df") <- nBs
  attr(BsData, "tgrid") <- tgrid
  attr(BsData, "Xgrid") <- Bs
  
  return(BsData)
}  
