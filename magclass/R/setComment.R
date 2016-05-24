
setComment <- function(object,nm=NULL) {
  getComment(object) <- nm
  return(object)
}