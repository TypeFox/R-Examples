is.transiogram <-
function(object) {
  if(!is(object, "transiogram")) return(FALSE)
  if(!prod(c("Tmat", "lags", "type") %in% names(object))) return(FALSE)
  if(!length(names(object)) %in% c(3, 4)) return(FALSE)
  if(!is.array(object$Tmat)) return(FALSE)
  if(diff(dim(object$Tmat)[1:2]) != 0) return(FALSE)
  if(!is.null(object$LOSE)) {
    if(!is.array(object$LOSE)) return(FALSE)
    if(diff(dim(object$LOSE)[1:2]) != 0) return(FALSE)
  }
  if(!is.numeric(object$lags)) return(FALSE)
  if(!is.character(object$type)) return(FALSE)
  if(object$type != "Empirical" && object$type != "Theoretical") return(FALSE)
  return(TRUE)
}

