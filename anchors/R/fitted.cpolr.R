#######################################################################
##
## Function: fitted.cpolr()
## Author  : Jonathan Wand <wand@stanford.edu>
##           http://wand.stanford.edu
## Created :  2007-02-01
##
#######################################################################

fitted.anchors.cpolr <- function( object, average = FALSE, unconditional = FALSE,  ...) {
  z <- fitted.cpolr( object$cpolr, object$rank, average, unconditional, ...)
  class(z) <- "fitted.anchors.cpolr"
  return(z)
}

fitted.cpolr <- function(object, anchors, average = FALSE, unconditional = FALSE,  ...) {

  ## only anchors anchorss make sense..
  if (!missing(anchors) && !(class(anchors) %in% c("anchors.rank","anchors.rank.type")))
    stop("Second argument, anchors, must be of class 'anchors.rank' or 'anchors.rank.type' \n")
  
  if (!missing(anchors) && class(anchors) == "anchors.rank" && class(anchors$rank) == "anchors.rank.type")
    anchors <- anchors$rank
  
  if (unconditional || missing(anchors)) {
    ## unconditional fitted values
    if (!missing(anchors)) {
      midx <- rownames(object$fitted.values) %in% rownames(anchors$span)
      mf2 <- object$fitted.values[midx,]
    } else {
      mf2 <- object$fitted.values
    }
  } else {
    ## conditional fitted values
    aidx <- rownames(anchors$span) %in% rownames(object$fitted.values)
    midx <- rownames(object$fitted.values) %in% rownames(anchors$span)
    if (sum(aidx) != sum(midx)) stop("fitted.cpolr: mismatch in rownames")

    mf <- object$fitted.values[midx,] *  matrix( as.numeric( anchors$weight[aidx,] > 0) , ncol=anchors$max)
    mf2 <- (mf / apply(mf,1,sum) )
    rownames(mf2) <- rownames(anchors$span)[aidx]
  }

  if (average)
    mf2 <-  apply(mf2,2,mean)
  
  class(mf2) <- "fitted.cpolr"
  return(mf2)
}
