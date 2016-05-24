fusedlasso.mod <- function(y, X, D, approx=FALSE,
                       maxsteps=2000, minlam=0, tol=1e-11,
                       ridge = T, eps=1e-8) {
  cl = match.call()

  if (missing(y)) stop("y is missing.")
  if (!is.numeric(y)) stop("y must be numeric.")
  if (length(y) == 0) stop("There must be at least one data point [must have length(y) > 1].")
  
  if (!is.null(X) && !is.matrix(X)) stop("X must be a matrix.")
  # Right now there is no check for X having full
  # column rank; this is for efficiency, we never
  # have to compute its pseudoinverse, otherwise
  
  
  # Make D a sparse matrix if it's not already
  #### 9 Oct
  #if (c(attributes(class(D))$package,"")[[1]] != "Matrix") {
  #  D = Matrix(D,sparse=TRUE)
  #}

  # Check dimensions
  if (is.null(X) && length(y)!=ncol(D)) stop("Dimensions don't match [length(y) != ncol(D)].")
  if (!is.null(X) && length(y)!=nrow(X)) stop("Dimensions don't match [length(y) != nrow(X)].")
  if (!is.null(X) && ncol(D)!=ncol(X)) stop("Dimensions don't match [ncol(X) != ncol(D)].")
  
  	p = dim(X)[2]
  	
    	if(ridge){
    		x = rbind(X, sqrt(eps) * c(rep(1,p)))
      		y2 = as.vector(c(y,0))
	      out = dualpathFusedX.mod(y2,x,D,approx,maxsteps,minlam,tol)
    	}
    	else out = dualpathFusedX.mod(y,X,D,approx,maxsteps,minlam,tol)


  out$call = cl
  
  return(out)
  
}
