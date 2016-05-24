#' Change a multivariate polynomial into a function.
#' 
#' Transforms an mpoly object into a function which can be
#' evaluated.
#' 
#' @param x an object of class mpoly
#' @param varorder the order in which the arguments of the function
#'   will be provided
#' @param vector whether the function should take a vector argument
#'   (TRUE) or a series of arguments (FALSE)
#' @param ... any additional arguments
#' @usage \method{as.function}{mpoly}(x, varorder = vars(x), vector
#'   = TRUE, ...)
#' @seealso \code{\link{plug}}
#' @export
#' @examples
#' 
#' p <- mp("x + 3 x y + z^2 x")
#' f <- as.function(p)
#' f(1:3) # -> 16
#' f(c(1,1,1)) # -> 5
#' 
#' f <- as.function(p, vector = FALSE)
#' f(1, 2, 3) # -> 16
#' f(1, 1, 1) # -> 5
#' 
#' f <- as.function(p, varorder = c('z','y','x'), vector = FALSE)
#' f(3, 2, 1) # -> 16
#' f(1, 1, 1) # -> 5
#' 
#' # for univariate mpolys, as.function() returns a vectorized function
#' # that can even apply to arrays
#' p <- mp("x^2")
#' f <- as.function(p)
#' f(1:10)
#' (mat <- matrix(1:4, 2))
#' f(mat)
#' 
#' 
#' p <- mp('1 2 3 4')
#' f <- as.function(p)
#' f(10) # -> 24
#' 
as.function.mpoly <- function(x, varorder = vars(x), vector = TRUE, ...){
	
  ## argument checking
  stopifnot(is.character(varorder))
  stopifnot(is.logical(vector))  	
  stopifnot(is.mpoly(x))
	
  if(!setequal(varorder, vars(x))){
    stop('varorder must contain all of the variables of x.',
      call. = FALSE)
  }
  
  ## determine the number of variables
  p <- length(vars(x))
 
  ## deal with constant polynomials
  if(is.constant(x)) return( function(.) unlist(x)[["coef"]] )
    
  ## univariate polynomial
  if(p == 1){
    mpoly_string <- suppressMessages(print.mpoly(x, stars = TRUE))
    mpoly_string <- chartr(vars(x), ".", mpoly_string)
    message("f(.) with . = ", vars(x))
    f <- function(){}
    formals(f) <- alist(. = )
    #expression(if(length(.) > 1) return(sapply(., f))),
    body(f) <- as.call(c(
      as.name("{"),
      expression(if(length(.) > 1){
        .[] <- sapply(., f)
        return(.)
      }),
      parse(text = mpoly_string)
    ))
    return(f)
  }
  
  ## general polynomials as a vector argument
  if(vector){
    mpoly_string <- suppressMessages(print.mpoly(x, stars = TRUE))
    mpoly_string <- paste(' ', mpoly_string, ' ', sep = '')
    for(k in 1:p){
      mpoly_string <- gsub(
        paste(' ', varorder[k], ' ', sep = ''),
        paste(' .[', k, '] ', sep = ''),
        mpoly_string
      )
      mpoly_string <- gsub(
        paste(' ', varorder[k], '\\*\\*', sep = ''),
        paste(' .[', k, ']**', sep = ''),
        mpoly_string
      )      
    }
    v <- paste('(', paste(varorder, collapse = ', '), ')', sep = '')
    message('f(.) with . = ', v)
    mpoly_string <- paste('function(.){', mpoly_string, '}')    
    return(eval(parse(text = mpoly_string)))
  }
  
  ## general polynomials as a bunch of arguments
  if(!vector){
    mpoly_string <- suppressMessages(print.mpoly(x, stars = TRUE))
    message('f(', paste(varorder, collapse = ', '), ')', sep = '')
    mpoly_string <- paste(
      'function(', 
      paste(varorder, collapse = ', '),
      '){', 
      mpoly_string, 
      '}',
      sep = ''
    )
    return(eval(parse(text = mpoly_string)))
  }  
  
}










as.function.bernstein <- function(x, ...){
 
  ## grab bernstein values
  k <- attr(x, "bernstein")$k
  n <- attr(x, "bernstein")$n
  
  ## return exp'd log function
  function(.) exp(lchoose(n, k) + k*log(.) + (n-k)*log(1-.))
  
}







as.function.jacobi <- function(x, ...){
  return(as.function.mpoly(x)) ## below is broken.
  
  ## grab bernstein values
  d <- attr(x, "jacobi")$degree
  k <- attr(x, "jacobi")$kind
  i <- attr(x, "jacobi")$indeterminate
  n <- attr(x, "jacobi")$normalized
  a <- attr(x, "jacobi")$alpha
  b <- attr(x, "jacobi")$beta
  
  ## return exp'd log function #
  #http://en.wikipedia.org/wiki/Jacobi_polynomials function(.)
  #pochhammer(a+1, d) / factorial(d) * hypergeo(-d, 1+a+b+d, a+1,
  #(1-.)/2)
  
}









