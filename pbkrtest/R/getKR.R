
getKR <- function (object, name = c("ndf", "ddf", "Fstat", "p.value", "F.scaling", "FstatU", "p.valueU", "aux")) 
{	
  stopifnot(is(object, "KRmodcomp"))
  if (missing(name) || is.null(name)){
    return(object$stats)
  } else {
    stopifnot(length(name <- as.character(name)) == 1)
    name <- match.arg(name)
    object$stats[[name]]
  }
}
