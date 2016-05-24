
#' Returns the constant part \code{b} of a linear (in)equality
#'
#' @example ../examples/editmatrixAttr.R
#' @export getb 
#' @seealso \code{\link{editmatrix}}
#'
#' @param E \code{\link{editmatrix}}
#' @return \code{numeric} vector \code{b}
getb <- function(E){
  if (!is.editmatrix(E)){
     stop("E has to be an editmatrix.")
  }
  E <- unclass(E)
  E[,ncol(E)]
}



#' Returns the coefficient matrix \code{A} of linear (in)equalities
#'
#' @example ../examples/editmatrixAttr.R
#' @export getA 
#' @seealso \code{\link{editmatrix}}
#' @aliases getA 
#'
#' @param E \code{\link{editmatrix}}
#' @return \code{numeric} matrix \code{A}
getA <- function(E){
  if ( is.editmatrix(E) ){
    unclass(E)[,-ncol(E),drop=FALSE]
  } else {
     stop("E has to be an editmatrix")
  }
}


#' Returns augmented matrix representation of edit set.
#'
#' For a system of linear (in)equations of the form \eqn{Ax \odot b}, \eqn{\odot\in\{<,\leq,=\}},
#' the matrix \eqn{A|b} is called the augmented matrix.
#'
#' @example ../examples/editmatrixAttr.R
#' @seealso \code{\link{editmatrix}} \code{\link{as.matrix.editmatrix}}
#'
#' @param E \code{\link{editmatrix}}
#' @return \code{numeric} matrix \code{A|b}
#' @export 
getAb <- function(E){
    if (!is.editmatrix(E))  stop("E has to be an editmatrix.")
    unclass(E)[,,drop=FALSE]
}

#' Returns the operator part of a linear (in)equality \code{editmatrix} E
#'
#' @export
#' @seealso \code{\link{editmatrix}}
#'
#' @example ../examples/editmatrixAttr.R
#' 
#' @param E \code{\link{editmatrix}}
#'
#' @return \code{character} vector with the (in)equality operators. 
getOps <- function(E){
  if (!is.editmatrix(E)){
     stop("E has to be an editmatrix.")
  }
  attr(E, "ops")
}


#' Check if an editmatrix is normalized
#'
#' @export
#' @seealso \code{\link{editmatrix}}
#'
#' @param E \code{\link{editmatrix}}
#'
#' @return TRUE when all comparison operators of \code{E} are in \{\code{<,<=,==}\}
isNormalized <- function(E){
  if (!is.editmatrix(E)){
     stop("Argument not of class editmatrix")
  }
  
  attr(E, "normalized") == TRUE ||
  all(getOps(E) %in% c("==","<","<="))
}

#' Normalizes an editmatrix
#'
#' An set of linear edits of the form \eqn{{\bf a}\cdot{\bf x}\odot b} with 
#' is called normalized when all  \eqn{\odot\in\{==,\leq,<\}} 
#'
#' @export
#' @seealso \code{\link{editmatrix}}
#'
#' @example ../examples/editmatrixAttr.R
#' 
#' @param E \code{\link{editmatrix}}
#'
#' @return If E was normalized, the original editmatrix is returned, otherwise 
#' a new normalized editmatrix will be returned 
normalize <- function(E){
  if (!is.editmatrix(E)) stop("Argument not of class editmatrix")
  if (isNormalized(E)){
     return(E)
  }
  
  A <- unclass(E)
  ops <- getOps(E)
  
  geq <- ops == ">="
  gt <- ops == ">"
  A[geq | gt,] <- -A[geq | gt,]
  ops[geq] <- "<="
  ops[gt] <- "<"      

  neweditmatrix(A, ops, normalized=TRUE)
}
