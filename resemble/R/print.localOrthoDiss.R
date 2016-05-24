#' @title Print method for an object of class \code{orthoDiss}
#' @description Prints the content of an object of class \code{orthoDiss}
#' @aliases print.localOrthoDiss
#' @usage \method{print}{localOrthoDiss}(x, ...)
#' @param x an object of class \code{localOrthoDiss} (returned by \code{orthoDiss} when it uses \code{local = TRUE}). 
#' @param ... arguments to be passed to methods (not yet functional).
#' @author Leonardo Ramirez-Lopez and Antoine Stevens
#' @export

######################################################################
# resemble
# Copyrigth (C) 2014 Leonardo Ramirez-Lopez and Antoine Stevens
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
######################################################################

## History:
## 09.03.2014 Leo     The tryMod function was removed

print.localOrthoDiss <- function(x,...){
  obj <- x
  if(is.list(obj))
  {
    obj$dissimilarity <- as.matrix(round(obj$dissimilarity, getOption("digits")))
    obj$dissimilarity[is.na(obj$dissimilarity)] <- "*"
    dm <- format(obj$dissimilarity, digits = getOption("digits"), justify = "right")
    print(list(n.components = object$n.components, loc.n.components = obj$loc.n.components, dissimilarity = noquote(dm)))
    cat("*: local non-neighbor sample")
  }
  if(is.matrix(obj))
  {
    object <- as.matrix(round(obj, getOption("digits")))
    object[is.na(obj)] <- "*"
    dm <- format(object, digits = getOption("digits"), justify = "right")
    print(dm, quote = FALSE)
    cat("*: local non-neighbor sample")
  }
}


"[.localOrthoDiss" <- function(x, rr, cl, drop = FALSE, ...){
  object <- x
  if(!is.logical(drop))
    drop <- FALSE
  class(object) <- NULL
  obj <- object[rr,cl, drop = drop]
  if(!drop)
    class(obj) <- "localOrthoDiss"
  return(obj)
}  
