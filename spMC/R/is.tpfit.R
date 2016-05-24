is.tpfit <-
function(object) {
  if(!is(object, "tpfit")) return(FALSE)
  if(!prod(c("coefficients", "prop", "tolerance") %in% names(object))) return(FALSE)
  if(length(names(object)) != 3) return(FALSE)
  if(!is.matrix(object$coefficients)) return(FALSE)
  if(!is.numeric(object$prop)) return(FALSE)
  if(!is.numeric(object$tolerance)) return(FALSE)
  return(TRUE)
}
