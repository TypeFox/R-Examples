#' Polynomial components
#' 
#' Compute quantities/expressions related to a multivariate 
#' polynomial.
#' 
#' @param x an object of class mpoly
#' @param ndx a subsetting index
#' @param varorder the order of the variables
#' @param order a total order used to order the terms
#' @param reduced if TRUE, don't include zero degrees
#' @param ... additional arguments
#' @return An object of class mpoly or mpolyList, depending on the
#'   context
#' @name components
#' @examples
#' (p <- mp("x y^2 + x (x+1) (x+2) x z + 3 x^10"))
#' p[2]
#' p[-2]
#' p[2:3]
#' 
#' LT(p)
#' LC(p)
#' LM(p)
#' 
#' multideg(p)
#' totaldeg(p)
#' monomials(p)
#' 
#' exponents(p)
#' exponents(p, reduce = TRUE)
#' lapply(exponents(p), is.integer)
#' 
#' homogeneous_components(p)
#' 





#' @rdname components
#' @export
`[.mpoly` <- function(x, ndx){
  x <- unclass(x)[ndx]
  class(x) <- "mpoly"
  x
}





#' @rdname components
#' @export
LT <- function(x, varorder = vars(x), order = "lex"){
  reorder.mpoly(x, varorder, order)[1]
}





#' @rdname components
#' @export
LC <- function(x, varorder = vars(x), order = "lex"){
  lt <- LT(x, varorder, order)
  unname(lt[[1]]["coef"])
}




#' @rdname components
#' @export
LM <- function(x, varorder = vars(x), order = "lex"){
  lt <- LT(x, varorder, order)
  lt[[1]]["coef"] <- 1
  lt
}





#' @rdname components
#' @export
multideg <- function(x, varorder = vars(x), order = "lex"){
  lt <- LT(x, varorder, order)
  coef_ndx <- which(names(lt[[1]]) == "coef")
  sparse_multideg <- lt[[1]][-coef_ndx]
  lt_vars <- names(sparse_multideg)
  nvars <- length(varorder)
  mdeg <- rep(0L, length = nvars)
  names(mdeg) <- varorder
  for (k in 1:length(lt_vars)) { 
    mdeg[lt_vars[k]] <- mdeg[lt_vars[k]] + sparse_multideg[k]
  }
  mdeg
}





#' @rdname components
#' @export
totaldeg <- function(x){
  if(!is.mpoly(x) && length(x) > 1){
    return(vapply(x, totaldeg, numeric(1)))
  }
  if(!is.mpoly(x)) stop("totaldeg requires an mpoly or mpolyList object.")
  max(vapply(exponents(x), sum, numeric(1)))
}





#' @rdname components
#' @export
monomials <- function(x){
  if(!is.mpoly(x)) stop("monomials requires an mpoly or mpolyList object.")
  xs <- lapply(1:length(x), function(ndx) `[.mpoly`(x,ndx))
  class(xs) <- "mpolyList"
  xs
}






#' @rdname components
#' @export
exponents <- function(x, reduced = FALSE){

  l <- lapply(x, function(term){
    fixed_term <- as.integer(term[-which(names(term) == "coef")])
    names(fixed_term) <- names(term[-which(names(term) == "coef")])
    fixed_term
  })  
  
  if(reduced) return(l)
  
  v <- vars(x)
  p <- length(v)
  tmp <- rep.int(0L, p)
  names(tmp) <- v
  
  lapply(l, function(exp){
    tmp <- tmp + exp[v]
    tmp[is.na(tmp)] <- 0
    tmp <- as.integer(tmp)
    names(tmp) <- v
    tmp
  })
}