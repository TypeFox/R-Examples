#######################################################################
##
## Function: fitted.anchors()
## Author  : Jonathan Wand <wand(at)stanford.edu>,
##           http://wand.stanford.edu
## Created :  2008-05-01
##
## MODIFIED:
##   
#######################################################################
fitted.anchors.rank <- function(object, ... ,
                                ties=c("omit","uniform","cpolr","minentropy"),
                                average = FALSE, unconditional=FALSE) {

  if (class(object) != "anchors.rank")
    stop("object must be of class anchors.rank")
  
  ties <- match.arg(ties)
  if (average) {

    if (ties=="omit") 
      out <- object$summary$scalar$Prop
    if (ties == "uniform") 
      out <- object$summary$uniform$Prop
    if (ties == "cpolr") 
      out <- fitted.cpolr(object$cpolr, object$rank, average=TRUE, unconditional=unconditional)     
    if (ties == "minentropy") 
      out <- summary.minimum.entropy(object$minentropy, average=TRUE)

  } else {

    if (ties=="omit")
      out <- object$rank$span[ object$rank$span[,1]==object$rank$span[,2] ,1]
    if (ties == "uniform") 
      out <- object$rank$weight
    if (ties == "cpolr") 
      out <- fitted(object$cpolr, object$rank, average=FALSE, unconditional=unconditional)     
    if (ties == "minentropy") 
      out <- summary(object$minentropy, average=FALSE)
    
  }

  return(out)
  
}

