outlierProbs <-
  ## Short form for generic function 
  function(object) UseMethod("outlierProbs")

outlierProbs.metaplus <- function(object) {
  if (!inherits(object, "metaplus"))
    stop("Use only with 'metaplus' objects.\n")
  if (object$random!="mixture")  
    stop("Use only with 'mixture' distribution.\n")
  outliers <- list(outlier.prob=object$outlier.prob,slab=object$slab)
  class(outliers) <- "outlierProbs"
  return(outliers)
}