#' Compute partial derivatives of a multivariate polynomial.
#' 
#' This is a deriv method for mpoly objects.  It does not call the
#' deriv function (from package stats).
#' 
#' @param expr an object of class mpoly
#' @param var character - the partial derivative desired
#' @param ... any additional arguments
#' @return An object of class mpoly or mpolyList.
#' @export
#' @examples
#' m <- mp('x y + y z + z^2')
#' deriv(m, 'x')
#' deriv(m, 'y')
#' deriv(m, 'z')
#' deriv(m, c('x','y','z'))
#' deriv(m, 'a')
#' is.mpoly(deriv(m, 'x'))
#' is.mpolyList( deriv(m, c('x','y','z')) )
deriv.mpoly <- function(expr, var, ...){

  if(missing(var)){
  	stop('var must be specified, see ?deriv.mpoly', call. = FALSE)
  }
	
  # argument checks	
  if(!is.mpoly(expr)){
    stop('expr must be of class mpoly.', call. = FALSE)
  }
  stopifnot(is.character(var))
  
  
  # if many vars provided
  if(length(var) > 1){
    mpolyList <- lapply(as.list(var), function(var){
      deriv(expr, var = var, ...)
    })
    class(mpolyList) <- 'mpolyList'  
    return(mpolyList)    
  }
  
	
  # if the variable is not in the polynomial, return zero
  if(!(var %in% vars(expr))){
    return( mpoly(list(c(coef = 0))) )
  }
	
  expr <- unclass(expr)
  
  # take derivative
  expr <- lapply(expr, function(v){
  	if(length(v) == 1) return(c(coef = 0))
    p <- length(v)
    if(!(var %in% names(v[1:p]))) return(c(coef = 0))
    v['coef'] <- unname(v[var]) * v['coef']
    v[var] <- v[var] - 1
    v
  })
  
  # return
  mpoly(expr)

}


#' Compute gradient of a multivariate polynomial.
#'
#' This is a wrapper for deriv.mpoly.
#'
#' @param mpoly an object of class mpoly
#' @seealso \code{\link{deriv.mpoly}}
#' @return An object of class mpoly or mpolyList.
#' @export
#' @examples
#' m <- mp('x y + y z + z^2')
#' gradient(m)
#' 
#' 
#' # gradient descent illustration using the symbolically
#' # computed gradient of the rosenbrock function
#' rosenbrock <- mp("(1 - x)^2 + 100 (y - x^2)^2")
#' fn <- as.function(rosenbrock)
#' (rosenbrock_gradient <- gradient(rosenbrock))
#' gr <- as.function(rosenbrock_gradient)
#' 
#' # visualize the function
#' library(ggplot2)
#' s <- seq(-5, 5, .05)
#' df <- expand.grid(x = s, y = s)
#' df$z <- apply(df, 1, fn)
#' ggplot(df, aes(x = x, y = y)) +
#'   geom_raster(aes(fill = z)) +
#'   scale_fill_continuous(trans = "log10")
#'   
#' # run the gradient descent algorithm using line-search
#' # step sizes computed with optimize()
#' current <- steps <- c(-3,-4)
#' change <- 1
#' tol <- 1e-5
#' while(change > tol){
#'   last  <- current
#'   delta <- optimize(
#'     function(delta) fn(current - delta*gr(current)),
#'     interval = c(1e-10, .1)
#'   )$minimum
#'   current <- current - delta*gr(current)
#'   steps   <- unname(rbind(steps, current))
#'   change  <- abs(fn(current) - fn(last))
#' }
#' steps <- as.data.frame(steps)
#' names(steps) <- c("x", "y")
#'   
#' # visualize steps, note the optim at c(1,1)
#' # the routine took 5748 steps
#' ggplot(df, aes(x = x, y = y)) +
#'   geom_raster(aes(fill = z)) +
#'   geom_path(data = steps, color = "red") +
#'   geom_point(data = steps, color = "red", size = .5) +
#'   scale_fill_continuous(trans = "log10")
#' 
#' # it gets to the general region of space quickly
#' # but once it gets on the right arc, it's terrible
#' # here's what the end game looks like
#' last_steps <- tail(steps, 100)
#' rngx <- range(last_steps$x); sx <- seq(rngx[1], rngx[2], length.out = 201)
#' rngy <- range(last_steps$y); sy <- seq(rngy[1], rngy[2], length.out = 201)
#' df <- expand.grid(x = sx, y = sy)
#' df$z <- apply(df, 1, fn) 
#' ggplot(df, aes(x = x, y = y)) +
#'   geom_raster(aes(fill = z)) +
#'   geom_path(data = last_steps, color = "red", size = .25) +
#'   geom_point(data = last_steps, color = "red", size = 1) +
#'   scale_fill_continuous(trans = "log10") 
#'   
#' 
gradient <- function(mpoly){
  deriv(mpoly, var = vars(mpoly))
}


