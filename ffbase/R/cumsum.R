#' Cumulative Sums, Products, and Extremes
#'
#' @method cumsum ff
#' @example ../examples/cumsum.R
#' @param x \code{ff} numeric vector or an object that can be coerced to one a numeric vector
#' @param ... other parameters passed on to chunk
#' @return An \code{ff} vector of the same length and type as x (after coercion), except that cumprod returns a numeric vector for integer input.\cr
#' An NA value in x causes the corresponding and following elements of the return value to be NA, as does integer overflow in cumsum (with a warning). 
#' @rdname cumsum.ff
#' @export
#' @export cumsum.ff
#' @seealso \code{\link{cumsum}}, \code{\link{cumprod}}, \code{\link{cummax}}, \code{\link{cummin}}
cumsum.ff <- function(x, ...){
  result <- ff::clone.ff(x, vmode = "double")
  
  i.last <- 0
  for (i in chunk(x, ...)){
    Log$chunk(i)
    cs <- cumsum(c(i.last, x[i]))
    i.last <- tail(cs,1)
    result[i] <- cs[-1]
  }
  result
}

#' @rdname cumsum.ff
#' @method cumprod ff
#' @export
cumprod.ff <- function(x, ...){
  result <- ff::clone.ff(x, vmode = "double")
  
  i.last <- 1
  for (i in chunk(x, ...)){
    Log$chunk(i)
    cs <- cumprod(c(i.last, x[i]))
    i.last <- tail(cs,1)
    result[i] <- cs[-1]
  }
  result
}

#' @rdname cumsum.ff
#' @method cummax ff
#' @export
cummax.ff <- function(x, ...){
  result <- ff::clone.ff(x, vmode = "double")
  
  i.last <- -Inf
  for (i in chunk(x, ...)){
    cs <- cummax(c(i.last, x[i]))
    i.last <- tail(cs,1)
    result[i] <- cs[-1]
  }
  result
}

#' @rdname cumsum.ff
#' @method cummin ff
#' @export
cummin.ff <- function(x, ...){
  result <- ff::clone.ff(x, vmode = "double")
  
  i.last <- Inf
  for (i in chunk(x, ...)){
    Log$chunk(i)
    cs <- cummin(c(i.last, x[i]))
    i.last <- tail(cs,1)
    result[i] <- cs[-1]
  }
  result
}
