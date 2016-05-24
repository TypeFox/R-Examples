#' Helper function to ensure that mpcross objects are in map order
#'
#' Orders markers within all components of mpcross object (founders, finals, etc.) to match the map order
#' @export
#' @param object Object of class \code{mpcross}
#' @return The original object, with all components reordered to match the map. Additional unmapped markers will appear after all mapped markers. 

maporder <- function(object)
{
  if (!inherits(object, "mpcross")) stop("Object must be of class mpcross")
  ## basically just make sure everything's in map order. 
  if (is.null(object$map)) {
	cat("No need to reorder, no map exists \n") 
	return(object)
  } 

  mapmrk <- unlist(lapply(object$map, names))
  allmrk <- c(mapmrk, setdiff(colnames(object$finals), mapmrk))
 
  m <- match(allmrk, colnames(object$finals))
  object$finals <- object$finals[,m]
  object$founders <- object$founders[,m]
  if (!is.null(object$rf)) {
  	m2 <- match(allmrk, colnames(object$rf$theta))
	object$rf$theta <- object$rf$theta[m2,m2]
	object$rf$lkhd <- object$rf$lkhd[m,m]
	object$rf$lod <- object$rf$lod[m,m]
  }
 
  if (inherits(object, "mpqtl")) 
	cat("Warning: May need to rerun QTL analysis because map has changed\n")
  if (inherits(object, "mpprob"))
	cat("Warning: May need to rerun probability computation because map has changed\n")

  object
}

