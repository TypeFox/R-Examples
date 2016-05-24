search.sphere <- function(k=9) {
  phi <- seq(0, 2*pi, length=k)
  theta <- seq(0, pi, length=k)
  
  angles <- data.frame(expand.grid(list(phi=phi, theta=theta)))
  angles <- subset(angles, (theta != 0) & (phi != 0))
  
  us <- with(angles, data.frame(U1=sin(theta)*cos(phi), U2=sin(theta)*sin(phi), U3=cos(theta)))
  us
}

sphereA <- function(A, theta, sphere) {
  ## gives a grid of rotation matrices with the same distance from rotation A
  R <- as.SO3(as.matrix(sphere), theta=rep(theta, length=nrow(sphere)))
  # multiplication isn't right
  X <- matrix(A, nrow=3)
  for (i in 1:nrow(R)) {
    R[i,] <- as.SO3(X %*% matrix(R[i,], nrow=3))   
  }
  R
}



L2.error <- function(sample, Shat) {
  sum(rot.dist(sample, Shat, method="intrinsic", p=2))
}

error.grid <- function(sample, Shat, theta=1, error, sphere) {
  rShat <- sphereA(Shat, theta, sphere)
  err <- vector(length=nrow(rShat))
  
  for (i in 1:nrow(rShat)) {
    R <- matrix(unlist(rShat[i, 1:9]), ncol=3)
    err[i] <- error(sample, R)
  }
  return(err)
}

#' Gradient optimization for rotation data
#' 
#' Gradient based optimization for user defined central orientation of a rotation sample.
#' 
#' @param sample sample of rotations.
#' @param error user defined function to observed distance between sample and estimate, has to have parameters for the sample and the estimate.
#' @param minerr minimal distance to consider for convergence.
#' @param start starting value for the estimation.
#' @param theta size of the grid considered.
#' @return list of 
#' \itemize{
#' \item \code{Shat} estimate of the main direction
#' \item \code{iter} number of iterations necessary for convergence
#' \item \code{theta} final size of the grid
#' \item \code{minerr} error used for convergence
#' \item \code{error} numeric value of total sample's distance from main direction
#' }
#' @export
#' @examples 
#' # minimize L1 norm:
#' L1.error <- function(sample, Shat) {
#'   sum(rot.dist(sample, Shat, method = "intrinsic", p = 1))
#' }
#' 
#' cayley.sample <- ruars(n = 10, rangle = rcayley, nu = 1, space = 'SO3')
#' SL1 <- gradient.search(cayley.sample, L1.error, start = id.SO3)
#' 
#' # visually no perceptible difference between median estimates from in-built function and 
#' # gradient based search (for almost all starting values)
#' 
#' \dontrun{
#' plot(cayley.sample, center=SL1$Shat, show_estimates="all")}

gradient.search <- function(sample, error, minerr =1e-5, start = mean(sample), theta=NULL) {
# 	if (length(start) == 1)
# 		Shat <- as.SO3(sample[start,])
# 	else if (all(dim(start) == 3))
# 		Shat <- as.SO3(start)

  Shat <- start
  sphere <- search.sphere()
  err <- error(sample, Shat)
	if (is.null(theta)) theta <- 0.5*err/nrow(sample)
  
	iter <- 0
	while (theta > minerr) {
	  rA <- sphereA(Shat, theta, sphere=sphere)
	  rA.error <- error.grid(sample, Shat, theta, error, sphere=sphere)
	  if (min(rA.error, na.rm=TRUE) < err) {
	    Shat <- as.SO3(matrix(unlist(rA[which.min(rA.error),]), ncol=3))
	  } else {
	    theta <- theta/2
	  }
	  err <- error(sample, Shat)
	  iter <- iter+1
	}

#	print(paste("iterations:",iter))
	return(list(Shat=Shat, iter=iter, theta=theta, minerr=minerr, error=err))
}

