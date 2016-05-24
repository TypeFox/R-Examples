lhsDesign <- function(n, dimension, randomized=TRUE, seed=NULL){
  
  # arguments : n         = number of points
  #             dimension = number of variables
  #  	    randomized = logical for randomized or centered points
  #             seed       = value of the random seed
  # output : a list containing the arguments and the LHS design
  
  
  # if no seed is provided in argument, choice of the seed for 'runif' and 'sample'
  if (is.null(seed)){
    seed <- as.numeric(Sys.time())
  }
  set.seed(seed)
  
  # Randomized LHS: U[0,1]-sampling of n x dim values
  if (randomized) ran = matrix(runif(n*dimension),nrow=n,ncol=dimension) 
  # Centered LHS
  else ran = matrix(0.5,nrow=n,ncol=dimension) 
    
  x = matrix(0,nrow=n,ncol=dimension)  # initializing matrix x
  
  for (i in 1:dimension) {
    idx = sample(1:n)        # vector of permutations of [1 to n]
    P = (idx-ran[,i]) / n    # vecteur of probabilities
    x[,i] <- P  }
  
  # Outputs:
  return(list(n=n,dimension=dimension,design=x,randomized=randomized,seed=seed))
}
