#######################################################################
##
## Function: summary.anchors.entropy
## Author  : Jonathan Wand <wand(at)stanford.edu>
##           http://wand.stanford.edu
## Created :  2007-02-02
##
## MODIFIED: 2008-05-01 : JW
## - flex sort options
#######################################################################
summary.anchors.combn <- function(object, ...,
                                    sort = c("max","estimated","minimum","interval","span"),
                                    digits=3) {

  sort <- match.arg( sort )

  idx <- switch(sort,
                max       = rev(order(object$Max)),
                estimated = rev(order(object$Est)),
                minimum   = rev(order(object$Min)),
                interval  = order(object$N.int),
                span      = order(object$Avg))
                
                
  
  cat("\nSummary of entropy and intervals by subsets of vignettes:\n\n")
  print( round( object[idx,], digits) )
}
