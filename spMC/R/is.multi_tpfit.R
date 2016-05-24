is.multi_tpfit <-
function(object) {
  if(!is(object, "multi_tpfit")) return(FALSE)
  if(!prod(c("coordsnames", "coefficients", "prop", "tolerance") %in% names(object))) return(FALSE)
  if(length(names(object)) != 5 - is.null(object$rotation)) return(FALSE)
  if(!is.list(object$coefficients)) return(FALSE)
  if(!prod(sapply(object$coefficients, is.list))) return(FALSE)
  if(!prod(sapply(object$coefficients, names) == "coefficients")) return(FALSE)  
  if(!prod(sapply(object$coefficients, function(xx) is.matrix(xx$coefficients))))  return(FALSE)
  if(!is.numeric(object$prop)) return(FALSE)
  if(!is.numeric(object$tolerance)) return(FALSE)
  if(!is.null(object$rotation)) {
    if(!is.numeric(object$rotation)) return(FALSE)
    if(length(object$rotation) != length(object$coefficients) - 1) return(FALSE)
  }
  return(TRUE)
}

