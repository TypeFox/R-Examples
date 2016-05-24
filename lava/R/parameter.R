##' @export
"parameter<-" <- function(x,...,value) UseMethod("parameter<-")

##' @export
"parameter<-.lvmfit" <- function(x,...,value) {
  parameter(Model(x),...) <- value
  return(x)
}


##' @export
"parameter<-.lvm" <- function(x,constrain,start,remove=FALSE,...,value) {
  if (inherits(value,"formula")) value <- all.vars(value)
  if (remove) {
      x$expar[value] <- NULL
      x$exfix[value] <- NULL
      x$attributes$parameter[value] <- NULL
      index(x) <- reindex(x)
      return(x)
      
  }
  if (!missing(start)) {
      if (length(start) != length(value)) stop("'start' and 'value' should be of the same lengths")
      start <- as.list(start)
      names(start) <- value
  } else {
      start <- as.list(rep(0,length(value))); names(start) <- value
  }
  if (!missing(constrain)) {
      newfix <- constrain
      if (!is.list(newfix)) newfix <- as.list(newfix)
  } else {
      newfix <- as.list(value);
  }
  names(newfix) <- value
  x$expar[value] <- start
  x$exfix[value] <- newfix
  index(x) <- reindex(x)
  x$attributes$parameter[value] <- TRUE
  return(x)
}

##' @export
parameter <- function(x,var,...) {
    if (missing(var)) return (names(unlist(x$attributes$parameter)))
    parameter(x,...) <- var
}
