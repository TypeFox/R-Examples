#' Bezier polynomials
#' 
#' Compute the Bezier polynomials of a given collection of points. 
#' Note that using \code{\link{as.function}} on the resulting Bezier
#' polynomials is made numerically stable by taking advantage of de
#' Casteljau's algorithm; it does not use the polynomial that is
#' printed to the screen.  See \code{\link{bezierFunction}} for
#' details.
#' 
#' @param ... either a sequence of points or a matrix/data frame of
#'   points, see examples
#' @param indeterminate the indeterminate of the resulting
#'   polynomial
#' @return a mpoly object
#' @author David Kahle
#' @seealso \code{\link{bezierFunction}}
#' @export
#' @examples
#' 
#' p1 <- c(0,  0)
#' p2 <- c(1,  1)
#' p3 <- c(2, -1)
#' p4 <- c(3,  0)
#' bezier(p1, p2, p3, p4)
#' 
#' 
#' points <- data.frame(x = 0:3, y = c(0,1,-1,0))
#' bezier(points)
#' 
#' 
#' points <- data.frame(x = 0:2, y = c(0,1,0))
#' bezier(points)
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' # visualize the bernstein polynomials
#' 
#' library(ggplot2); theme_set(theme_bw())
#' 
#' s <- seq(0, 1, length.out = 101) 
#' 
#' ## example 1
#' points <- data.frame(x = 0:3, y = c(0,1,-1,0))
#' (bezPolys <- bezier(points))
#' 
#' f <- as.function(bezPolys)
#' df <- as.data.frame(f(s))
#' 
#' ggplot(aes(x = x, y = y), data = df) + 
#'   geom_point(data = points, color = "red", size = 8) +
#'   geom_path(data = points, color = "red") +
#'   geom_path()
#'   
#'   
#'   
#'   
#' ## example 1 with weights
#' f <- as.function(bezPolys, weights = c(1,5,5,1))
#' df <- as.data.frame(f(s))
#' 
#' ggplot(aes(x = x, y = y), data = df) + 
#'   geom_point(data = points, color = "red", size = 8) +
#'   geom_path(data = points, color = "red") +
#'   geom_path()
#'   
#'   
#'   
#'   
#'   
#' ## example 2
#' points <- data.frame(x = 0:2, y = c(0,1,0))
#' (bezPolys <- bezier(points))
#' f <- as.function(bezPolys)
#' df <- as.data.frame(f(s))
#' 
#' ggplot(aes(x = x, y = y), data = df) + 
#'   geom_point(data = points, color = "red", size = 8) +
#'   geom_path(data = points, color = "red") +
#'   geom_path()
#' 
#' 
#' 
#' 
#' ## example 3
#' points <- data.frame(x = c(-1,-2,2,1), y = c(0,1,1,0))
#' (bezPolys <- bezier(points))
#' f <- as.function(bezPolys)
#' df <- as.data.frame(f(s))
#' 
#' ggplot(aes(x = x, y = y), data = df) + 
#'   geom_point(data = points, color = "red", size = 8) +
#'   geom_path(data = points, color = "red") +
#'   geom_path()
#'   
#'   
#'   
#'   
#' ## example 4
#' points <- data.frame(x = c(-1,2,-2,1), y = c(0,1,1,0))
#' (bezPolys <- bezier(points))
#' f <- as.function(bezPolys)
#' df <- as.data.frame(f(s))
#' 
#' ggplot(aes(x = x, y = y), data = df) + 
#'   geom_point(data = points, color = "red", size = 8) +
#'   geom_path(data = points, color = "red") +
#'   geom_path()
#'   
#'   
#'   
#'   
#' ## example 5
#' qplot(speed, dist, data = cars)
#' 
#' s <- seq(0, 1, length.out = 201) 
#' p <- bezier(cars)
#' f <- as.function(p)
#' df <- as.data.frame(f(s))
#' qplot(speed, dist, data = cars) +
#'   geom_path(data = df, color = "red")
#' 
#' # the curve is not invariant to permutations of the points
#' # but it always goes through the first and last points
#' permute_rows <- function(df) df[sample(nrow(df)),]  
#' p <- bezier(permute_rows(cars))
#' f <- as.function(p)
#' df <- as.data.frame(f(s))
#' qplot(speed, dist, data = cars) +
#'   geom_path(data = df, color = "red")
#' 
#' 
#' 
bezier <- function(..., indeterminate = "t"){  
  
  ## grab input
  dots <- as.list(match.call(expand.dots = FALSE))$"..."
  names(dots) <- as.character(match.call(expand.dots = FALSE)$"...")
  dots <- lapply(dots, eval)
  
  ## parse input into a data frame
  if(length(dots) == 1 && (is.data.frame(dots[[1]]) || is.matrix(dots[[1]]))){    
    points <- dots[[1]]
    if(is.matrix(points)) points <- as.data.frame(points)
    if(is.null(names(points))) names(points) <- paste0("x", 1:ncol(points))
  } else if(
    all(vapply(dots, is.numeric, logical(1))) &&
    all(vapply(dots, length, numeric(1)) == length(dots[[1]]))
  ){
    points <- as.data.frame(do.call(rbind, dots))    
    if(is.null(names(dots[[1]]))){
      names(points) <- paste0("x", 1:ncol(points))
    } else {
      names(points) <- names(dots[[1]])
    }
    row.names(points) <- NULL
  }
  

  ## make polynomial
  n <- nrow(points) 
  bernPolys <- bernstein(0:(n-1), n-1, indeterminate)
  
  
  ## initialize bezPoly
  d <- ncol(points)
  bezPoly <- mp(rep("0", d))
  
  for(k in 1:n){
    ## convert vector of numerics into mpolyList
    polyWeights <- lapply(unname(points[k,]), function(.) mpoly(list(c(coef = .))))
    class(polyWeights) <- "mpolyList"
    
    ## duplicate the bern basis function to the same length
    basis <- replicate(d, bernPolys[[k]], simplify = FALSE)
    class(basis) <- "mpolyList"
    
    ## accumulate onto bezPoly
    bezPoly <- bezPoly + polyWeights * basis
  }
  
  ## class with meta-data
  class(bezPoly) <- c("bezier", "mpolyList")
  attr(bezPoly, "bezier") <- list(points = points)
  
  ## return
  bezPoly
}




































#' Bezier function
#' 
#' Compute the Bezier function of a collection of polynomials.  By
#' Bezier function we mean the Bezier curve function, a parametric
#' map running from t = 0, the first point, to t = 1, the last
#' point, where the coordinate mappings are linear combinations of
#' Bernstein polynomials.
#' 
#' The function returned is vectorized and evaluates the Bezier
#' curve in a numerically stable way with de Castlejau's algorithm
#' (implemented in R).
#' 
#' @param points a matrix or data frame of numerics.  the rows
#'   represent points.
#' @param weights the weights in a weighted Bezier curve
#' @return function of a single parameter
#' @author David Kahle
#' @references \url{http://en.wikipedia.org/wiki/Bezier_curve}, 
#'   \url{http://en.wikipedia.org/wiki/De_Casteljau's_algorithm}
#' @seealso \code{\link{bezier}}
#' @export
#' @examples
#' 
#' library(ggplot2); theme_set(theme_bw())
#' 
#' 
#' t <- seq(0, 1, length.out = 201)
#' points <- data.frame(x = 0:3, y = c(0,1,-1,0))
#' 
#' 
#' f <- bezierFunction(points)
#' df <- as.data.frame(f(t))
#' 
#' ggplot(aes(x = x, y = y), data = df) + 
#'   geom_point(data = points, color = "red", size = 8) +
#'   geom_path(data = points, color = "red") +
#'   geom_path()
#'   
#'   
#'   
#'   
#' f <- bezierFunction(points, weights = c(1,5,5,1))
#' df <- as.data.frame(f(t))
#' 
#' ggplot(aes(x = x, y = y), data = df) + 
#'   geom_point(data = points, color = "red", size = 8) +
#'   geom_path(data = points, color = "red") +
#'   geom_path()
#'   
#'   
#'   
#'   
#' f <- bezierFunction(points, weights = c(1,10,10,1))
#' df <- as.data.frame(f(t))
#' 
#' ggplot(aes(x = x, y = y), data = df) + 
#'   geom_point(data = points, color = "red", size = 8) +
#'   geom_path(data = points, color = "red") +
#'   geom_path()
#'   
#'   
#'   
#'   
#'   
#'   
#'   
#'   
#'   
#'   
#' 
bezierFunction <- function(points, weights = rep(1L, nrow(points))){
  
  n <- nrow(points)
  degree <- n-1
  points <- cbind(1, as.matrix(points))
  points <- weights * points
  
  combineTwo <- function(vec, t) t*vec[-length(vec)] + (1-t)*vec[-1]  
  
  singlePointBezier <- function(.){
    for(i in 1:degree) points <- apply(points, 2, combineTwo, t = .)
    points[-1] / points[1]
  }  
  
  function(.){
    if(length(.) > 1) return(t(sapply(., singlePointBezier)))
    singlePointBezier(.)
  }
}










#' @export 
#' @rdname as.function.mpolyList
as.function.bezier <- function(x, ...) bezierFunction(attr(x, "bezier")$points, ...)








