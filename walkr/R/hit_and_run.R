#' Hit and Run 
#' 
#' This function provides a wrapper 
#' for the har function of the 
#' hit-and-run package
#' 
#'@param A is the lhs of Ax <= b
#'@param b is the rhs of Ax <= b
#'@param x0 is the starting point (a list of points)
#'@param points is the number of points we want to sample
#'@param thin every thin-th point is stored
#'@param burn the first burn points are deleted
#'@param chains is the number of chains we run
#' 
#' @return a list of chains of the sampled points, each chain
#'         being a matrix object with each column as a point
#' 
#' @examples
#' \dontrun{
#' ## note that this Ax <= b is different from Ax=b that the 
#' ## user specifies for walkr (see transformation section in vignette) 
#' hit_and_run(A = A, b = b, x0, points = 100, thin = 1, burn = 0, chains = 1) 
#' } 
#' @importFrom hitandrun har 


hit_and_run <- function(A, 
                        b, 
                        x0, 
                        points, 
                        thin = 1,
                        burn = 0,
                        chains = 1) {
  
  ## thin needs to be a multiple of the total number of points 
  
  stopifnot(points %% thin == 0)
  
  
  ## initialize the list of constraints
  ## should also do some checking
  
  constr <- list(constr = A, rhs = b, dir = rep("<=", nrow(A)))
  
  ## getting samples
  
  answer <- list()
  
  for (j in 1:chains) {
    
    ## total points is : points * thin * 1/(1-burn) / chains 
    ## because burn-in is a percentage, we must take the CEILING function
    ## to sample more than we need (in the case where dividing by 1-burn
    ## does not return an integer)
    
    total.points <- ceiling( (points / chains)  * thin * (1/(1-burn))) 
    
    result <- t(hitandrun::har(x0[[j]], constr, N = total.points, 
                               thin = 1)$samples)
    
    ## NEED TO HANDLE THE CASE WHEN ALPHA IS JUST 1 DIMENSIONAL
    
    if(dim(result)[1] == 1) {
      
      ## first, delete out the number of points that we want to burn
      ## second, only take every thin-th point
      
      ## we take the floor function because we took the ceiling above
      ## so in the case that multiplying by burn doesn't result in an integer
      ## we returning the correct number of points
      
      result <- matrix(result[, (floor(burn*total.points)+1) : total.points], nrow = 1)
      result <- matrix(result[ , (1:(points/chains))*thin], nrow = 1)
      
    }
    
    else {
      
      ## first, delete out the number of points that we want to burn
      ## second, only take every thin-th point
      
      ## same as above
      
      result <- result[, (floor(burn*total.points)+1) : total.points]
      result <- result[ , (1:(points/chains))*thin]
      
    }
    
    answer[[j]] <- result
  }
    
  return(answer)
  
} 