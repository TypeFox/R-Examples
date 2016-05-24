##***********************************************************************
## methods to compute the spacings for diffrent objects.
##
## First sort in descending order, then diff and possibly weight.
##
## Author: Yves Deville <deville.yves@alpestat.com> 
##
##***********************************************************************

spacings <- function(object, ...) {
  UseMethod("spacings")
}

spacings.numeric <- function(object, wExp = TRUE, ...) {
  if(is.null(object)) return(NULL)  
  if (length(object) > 1L) {
    sp <- -diff(sort(object, decreasing = TRUE))
    if (wExp) sp <- sp * (1:length(sp))
  } else sp <- numeric(0)
  sp
}

##=======================================================================
## to be used with a data.frame with a "block" variable
## such as MAXdata of OTSdata
##=======================================================================

spacings.data.frame <- function(object, varName, wExp = TRUE, ...) {
  ff <- function(x) {
    if (length(x) > 1) return(-diff(sort(x, decreasing = TRUE)))
    else return(numeric(0))
  }
  ##
  sp <- tapply(X = object[ , varName], INDEX = object[ , "block"], FUN = ff)
  ## weight the spacings
  if (wExp) sp <- lapply(sp, function(x) { x * (1L:length(x)) })
  sigmaHat <- mean(unlist(sp))
  ## attr(sp, "weights") 
  sp
}

##=======================================================================
## the data.frame in a "Rendata" object can be choosed with 'type'
##=======================================================================

spacings.Rendata <- function(object,
                             type = c("MAX", "OTS", "OT"),
                             wExp = TRUE, ...) {
  
  type <- match.arg(type)
  varName <- object$info$varName
  
  if (type != "OT") {
    res <- spacings(object[[sprintf("%sdata", type)]], varName = varName,
                    wExp = wExp)
  } else {
    res <- spacings.numeric(object = object[["OTdata"]][ , varName], wExp = wExp)
  }
  
  res
  
}


if (FALSE) {
  rd <- Garonne
  res <- spacings(object = rd$MAXdata, varName = rd$info$varName)
  res1 <- spacings(rd, type = "OT")
  res2 <- spacings(rd, type = "MAX")
}

