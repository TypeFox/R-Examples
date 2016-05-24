#' start point
#' 
#' Given \eqn{Ax \le b}, which defines a convex
#' polytope, this function picks n random
#' starting "center" points using linear programming. 
#' 
#' @param A is the lhs of \eqn{Ax \le b}
#' @param b is the rhs of \eqn{Ax \le b}
#' @param n is the number of points we want to return
#' @param average is the number of boundary points we want to take the average
#'   of
#' 
#' @return a matrix, with each column as a point
#' 
#' @importFrom limSolve lsei
#' 
#' @examples
#' \dontrun{
#' ## note that this Ax <= b is different from Ax=b that the 
#' ## user specifies for walkr (see transformation section in vignette) 
#' start_point(A = A, b = b, n = 1, average = 10) 
#' }


start_point <- function(A, 
                        b, 
                        n = 1, 
                        average = 10) {
  
  ## initialize the return matrix, each column is 1 point
  ## the points have dimension equal to the rows of A
  
  result <- matrix(ncol = n, nrow = ncol(A))
 
  ## creating every point
  
  for (i in 1:n) {
  
    ## initialize local variable to store boundary points
    ## so that we could take an average in the end
    
    new_x0 <- numeric()
  
    ## average is the number of boundary points we want to take the average of 
    ## the higher this number is, the more likely our n points will be closer to
    ## each other
    
    ## the basic idea for having "average" number of points is that the linear
    ## programming method we use below in general finds points that randomly lie
    ## at the bounds of our convex polytope. And thus, if the take "average"
    ## number of points, and then take the average of the boundary points, we
    ## can get a good approximation of the "analytic center" of the space 
    ## however, I am unsure of how many points to average is good
    
    ## it is true that the more points we take, the closer it will be to the 
    ## true mean of whatever distribution this linear programming algorithm 
    ## draws from however, if we want a diversified sample of starting points,
    ## it's probably a good idea to have the "average" on the lower side, so
    ## that they are not close to each other. Also, this number should increase
    ## as the dimension of the problem increases (i.e. number of rows in A)
    
    for(j in 1:average) {
        
        ## these two lines randomize the objective function
        
        objfunc <- matrix(sample(stats::runif(1,-1, 1), ncol(A), replace = TRUE),
                          nrow = 1, ncol = ncol(A))
        const <- stats::runif(1, -1, 1)
        
        ## suppresswarnings because we don't specific a equality constraint,
        ## we only care about Ax <= b (lsei gives a warning, but it's fine)
        
        ## in terms of the lsei, it actually takes in Gx >= h
        ## thus, we pass in -Ax >= -b , which is equivalent to Ax <= b
        
        new_x0 <- rbind(new_x0, tryCatch (    
          suppressWarnings(
            
            ## + 0.00001 to avoid the optimization not starting
            
            limSolve::lsei(A = objfunc, B = const, G = -1* A, H = -1 * b + 0.000001)[[1]]),
          error = function(c) stop("The inequality constraints cannot be all satisfied!
                                   Sampling in this solution space is not possible!")
      ))
    }
    
    ## take the average of those boundary points ==> to obtain the center
    
    result[, i] <- apply(new_x0,  2, mean)
 
  }
  
  ## if user only wants 1 point
  ## then we should return the result as a vector for 
  ## the user's convenience
  
  if(n == 1) {
    return(as.vector(result))
  }
  
  else {
    return(result)
  }
  
}