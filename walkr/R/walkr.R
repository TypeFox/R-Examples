#' The Walkr Function 
#' 
#' Given \eqn{Ax = b}, \code{walkr} samples points from the intersection of 
#' \eqn{Ax = b} with the n-simplex (\eqn{\sum x = 1}, \eqn{x_i \ge 0}). The 
#' Ax = b must be underdetermined, otherwise there is an unique solution
#' and there will be no sampling. 
#' 
#' Before the sampling, walkr internally performs the affine transformation
#' which takes the complete solution of \eqn{Ax = b} and that intersected
#' with the unit simplex into a space parametrized by coefficients, which 
#' we call the "alpha-space". The specific set of procedures taken is 
#' written in detail in the vignette. Essentially, the space is transformed,
#' the sampling takes place in the transformed space, and in the end 
#' \code{walkr} transforms back into the original coordinate system 
#' and returns the result. This transformation is affine, so the uniformity 
#' and mixing properties of the MCMC algorithms are not affected. 
#' The current MCMC sampling methods supported are "hit-and-run", 
#' "dikin" and "optimized-dikin" (a Rcpp boosted version for speed). 
#' 
#' 1) Hit-and-run is computationally less expensive and also 
#' guarantees uniformity asympotically with complexity of O(n^3)
#' points with respect to dimension n. However, in real practice, 
#' as dimensions ramp up, the mixing of hit-and-run is poor compared 
#' to Dikin. Thus, a lot of thinning would be needed as dimension
#' ramps up.
#' 
#' 2) Dikin Walk is a nearly uniform method known for its very strong 
#' mixing properties. However, each Dikin step is much more 
#' computationally expensive than hit-and-run, so it takes more time 
#' to sample every point. Thus, the "dikin" method
#' uses RcppEigen to speed up the core computationally 
#' expensive operations in the algorithm.
#' 
#' @param A is the lhs of the matrix equation A
#' @param b is the rhs of the matrix equation b
#' @param points is the number of points we want to sample
#' @param method is the MCMC sampling method. Please enter "hit-and-run", "dikin", or "
#'        optimized-dikin"
#' @param thin every thin-th point is stored
#' @param burn the first burn points are deleted
#' @param chains is the number of chains we run
#' @param ret.format is the format in which walkr returns the answer. Please
#'        enter "list" (of chains) or "matrix". 
#' 
#' @return Either a list of chains (with each chain as a matrix of points) or 
#'         a matrix containing all the points. Each column is a point sampled. 
#'   
#' @examples
#' ## 4D constraint
#' A <- matrix(c(2,0,1,3), ncol = 4)
#' b <- 0.5
#' sampled_points <- walkr(A = A, b = b, points = 100, method = "dikin")                
#'                                                
#' @export 
#' 

walkr <- function(A, 
                  b, 
                  points, 
                  method = "dikin",
                  thin = 1,
                  burn = 0.5,
                  chains = 1,
                  ret.format = "matrix") {
  
  ## 0. Doing some checking here
  if(!is.matrix(A)) {
    stop("A needs to be a matrix")
  }
  
  if(!(ncol(A) > nrow(A))) {
    stop("A must be underdetermined")
  }
  
  if(nrow(A) != length(b)) {
    stop("Dimensions of A and b don't match")
  }
  
  if(!is.numeric(b)) {
    stop("b needs to be a numeric vector")
  }
  
  if(!is.numeric(A)) {
    stop("A needs to contain numbers only")
  }
  
  if(! (method %in% c("hit-and-run", "dikin"))) {
    stop("Method must be hit-and-run or dikin")
  } 
  
  if( (points < 1) | (thin < 1) | (chains < 1) ) {
    stop("points, thin, chains must be geq 1")
  }
  if(burn < 0) {
    stop("burn must be non-negative")
  }
  
  if(burn >= 1){
    stop("burn must be smaller than 100 percent")
  }
  
  if(points %% chains != 0){
    stop("chains must be a multiple of points")
  }
  
  ## the positive ( > 0) part was already checked above there
  
  if(points %% 1 != 0){
    stop("points must be a positive integer")
  }
  
  if(thin %% 1 != 0){
    stop("thin must be a positive integer")
  }
  
  if(chains %% 1 != 0){
    stop("chains must be a positive integer")
  }
  if(!(ret.format %in% c("list", "matrix"))){
    stop("ret.format must be either list or matrix")
  }
  
  ## augment the simplex constraints 
  
  aug_A <- rbind(A, matrix(1, ncol = ncol(A), nrow = 1))
  aug_b <- c(b, 1)

  ## 1. regardless of method, we need to perform the affine transformation which
  ## takes us from x-space (Ax = b) into the alpha-space in which the polytope
  ## described is Ax <= b (denoted below as new_A, new_b)
  ## From there, we could perform the sampling
  
  ## break it up into null space(homogeneous) and particular solution
  
  ## the user enters Ax = b (and we assume they want to intersect the solution
  ## space of that with the simplex) in order for hit-and-run / dikin to sample,
  ## we need a convex polytope. Therefore, we must perform the affine transformation (described in the vignette) which 
  ## brings us from x-space into alpha-space. Specifically, we want to 
  ## 1. represent the solution space in terms of basis 
  ## 2. rewrite it in the generic form of Ax <= b 
  
  ## First, we tag on the simplex equality constraints as an extra row in A.
  ## Next, we know the complete solution is written as v_particular +
  ## homogeneous, where homogeneous is an infinite set parameterized by alpha.
  ## An important note here is that the homogeneous solution returned by
  ## MASS::Null is an orthonormal basis. Since we have this alpha
  ## parametrization, we can write it as v_p + homogeneous %*% alpha >= 0 Then:
  
  ## Homogeneous %*% alpha >= -vp
  ## -Homogeneous %*% alpha <= vp    (Ax <= b form)
  
  z <- complete_solution(A = aug_A, b = aug_b)
  
  ## need the particular and homogeneous because in the end 
  ## we want to transform back in to "x-space"
  
  particular  <- z$particular
  homogeneous <- z$homogeneous
  
  
  new_A <- -homogeneous
  new_b <- particular
  
  
  ## 2. Find starting point within convex polytope
  x0 <- list()
  
  for (q in 1:chains) {
    x0[[q]] <- start_point(A = new_A, b = new_b, n = 1, average = 20)    
  }
  
  ## 3. The sampling
  
  if(method == "dikin") {
    
    ## sampling in alpha space
    ## n = n - 1 because dikin takes starting point as the 1st sampled point
      
    alphas <- dikin_walk(A = new_A, b = new_b, points = points, r = 1, 
                         x0 = x0, thin = thin, burn = burn, chains = chains)
    
    ## <convert back into x-space> there are two lambda functions here 1)
    ## mapping is the function to be applied to the individual columns of each
    ## of the chains (matricies). It acts of a matrix of alphas. mapping is then
    ## lapplied to each of the element in the list 2) the homogeneous x +
    ## particular function is acting on each column of the matrix, so that we
    ## can convert the alphas back into x-space
    
    
    mapping <- function(alpha_matrix) {apply(alpha_matrix, 2, function(x) { homogeneous %*% x + particular  })}   
    answer <- lapply(alphas, mapping)
    
    ## here, we calculate the rhats from the individual chains 
    ## we warn the user if any of the rhats for the parameters is above 1.1
    
    rhats <- calc_rhat(answer)
    rhats <- ifelse(!is.na(rhats), rhats, 1)
    
    if( any(rhats > 1.1) ) {
      warning("there are parameters with rhat > 1.1, you may want to run your chains for longer")
    }
  }
  
  else if (method == "hit-and-run") {
  
    
    ## to see what's going on here, see the comments above for dikin
    ## it is doing the same thing
    
    alphas <- hit_and_run(A = new_A, b = new_b, x0 = x0, points = points, 
                            thin = thin, burn = burn, chains = chains)
    
    mapping <- function(alpha_matrix) {apply(alpha_matrix, 2, function(x) { homogeneous %*% x + particular  })}
    
    answer <- lapply(alphas, mapping)
    rhats <- calc_rhat(answer)
    
    rhats <- calc_rhat(answer)
    rhats <- ifelse(!is.na(rhats), rhats, 1)
    
    if( any(rhats > 1.1) ) {
      warning("there are parameters with rhat > 1.1, you may want to run your chains for longer")
    }
  
  }
  
  else{
    stop("Sampling method must be \"hitandrun\" or \"dikin\".")
  }
  
  if(ret.format == "matrix") {
    ## after we perform all the checking above we return the chains as a matrix,
    ## with the sampled points as a row thus, here we column bind each of the
    ## chains, and then return the final result
    
    mat_answer <- answer[[1]]
    
    ## for loop breaks if chains == 1, so we need the if to check it
    
    if(chains > 1) {
      for(i in 2:chains) {
        
        mat_answer <- cbind(mat_answer, answer[[i]])
        
      }    
    }
  
  
    return(mat_answer)
  }
  
  ## since checking for correct input was done above, we
  ## can assume ret.format is 1 of the 2 correct ones
  
  else {
    
    ## list of chains
    
    return(answer)
  }

}