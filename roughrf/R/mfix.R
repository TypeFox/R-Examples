mfix <- function(x, mmmm) {
  
  m4m = function(x) {
    c(median(x, na.rm = T), mean(x, na.rm = T), min(x, na.rm = T), 
      max(x, na.rm = T))[mmmm]
  }
  
  mmfix <- function(object, ...) UseMethod("mmfix")
  
  mmfix.data.frame <- function(object, ...) {
    isfac <- sapply(object, is.factor)
    isnum <- sapply(object, is.numeric)
    if (any(!(isfac | isnum))) 
      stop("mfix only works for numeric or factor")
    roughfix <- function(x) {
      if (any(is.na(x))) {
        if (is.factor(x)) {
          freq <- table(x)
          x[is.na(x)] <- names(freq)[which.max(freq)]
        } else {
          x[is.na(x)] <- m4m(x)
        }
      }
      x
    }
    object[] <- lapply(object, roughfix)
    object
  }
  
  mmfix.default <- function(object, ...) {
    if (!is.atomic(object)) 
      return(object)
    d <- dim(object)
    if (length(d) > 2) 
      stop("can't handle objects with more than two dimensions")
    if (all(!is.na(object))) 
      return(object)
    if (!is.numeric(object)) 
      stop("mfix can only deal with numeric data.")
    if (length(d) == 2) {
      hasNA <- which(apply(object, 2, function(x) any(is.na(x))))
      for (j in hasNA) object[is.na(object[, j]), j] <- m4m(object[,j])
    } else {
      object[is.na(object)] <- m4m(object)
    }
    object
  }
  mmfix(x)
}
