#' Swap polynomial indeterminates
#' 
#' Swap polynomial indeterminates
#' 
#' @param p polynomial
#' @param variable the variable in the polynomial to replace
#' @param replacement the replacement variable
#' @return a mpoly object
#' @author David Kahle
#' @export
#' @examples
#' 
#' (p <- mp("(x + y)^2"))
#' swap(p, "x", "t")
#' 
#' ## the meta data is retained
#' (p <- bernstein(3, 5))
#' (p2 <- swap(p, "x", "t"))
#' is.bernstein(p2)
#' 
#' (p <- chebyshev(3))
#' (p2 <- swap(p, "x", "t"))
#' is.chebyshev(p2)
#' 
#' 
swap <- function(p, variable, replacement){


  ## arg check
  stopifnot(is.mpoly(p))

  ## determine variables
  vars <- vars(p)

  ## if constant, return
  if(is.constant(p)) return(p)
  

  ## don't allow replacement of one variable in the presence of many when
  ## the variable is already a part of the polynomial, since
  ## mpoly won't be run again to determine if a variable is duplicated
  stopifnot(variable %in% vars)
  if(replacement %in% vars && length(vars) > 1){
    stop("the replacement value cannot be a variable in the polynomial, try plug.", call. = FALSE)
  }
  
  
  ## swapping in polynomial
  pOut <- unclass(p)
  pOut <- lapply(p, function(v){
    names(v)[names(v) == variable] <- replacement
    v
  })
  class(pOut) <- "mpoly"
  
  ## custom for special polynomials
  if(is.bernstein(p)){
    class(pOut) <- c("bernstein", "mpoly")
    attr(pOut, "bernstein") <- attr(p, "bernstein")
    attr(pOut, "bernstein")["indeterminate"] <- replacement
  }
  
  if(is.chebyshev(p)){
    class(pOut) <- c("chebyshev", "mpoly")
    attr(pOut, "chebyshev") <- attr(p, "chebyshev")
    attr(pOut, "chebyshev")["indeterminate"] <- replacement
  }
  
  ## return
  pOut
}



