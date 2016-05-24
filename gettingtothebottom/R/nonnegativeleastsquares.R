#' Generate random nonnegative mixture components
#' 
#' \code{generate_nnm} Function to random nonnegative mixture components
#' 
#' @param n Number of samples
#' @param p Number of components
#' @param seed Random seed
#' 
#' @export
#' 
#' @examples
#' n <- 1e3
#' p <- 10
#' nnm <- generate_nnm(n,p)
#' 
generate_nnm <- function(n,p,seed=12345) {
  set.seed(seed)
  
  t <- seq(-10,10,length.out=n)
  # Choose some centers
  mu <- double(p)
  mu <- seq(-7,7,length.out=p)
  mu <- mu + rnorm(p)
  
  # Choose some widths
  sigma <- rgamma(p,shape=0.75)
  X <- matrix(NA,n,p)
  for (i in 1:p) {  
    X[,i] <- dnorm(t,mean=mu[i],sd=sigma[i])
  }
  return(list(X=X,t=t,mu=mu,sd=sigma))
}

#' Nonnegative Least Squares via MM
#' 
#' \code{nnls_mm} Iteratively computes the solution to the nonnegative least squares problem via a majorization-minimization algorithm.
#' 
#' @param y Nonnegative response
#' @param X Nonnegative design matrix
#' @param b Nonnegative initial regression vector
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance for convergence
#' 
#' @export
#' 
#' @examples
#' set.seed(12345)
#' n <- 100
#' p <- 3
#' X <- matrix(rexp(n*p,rate=1),n,p)
#' b <- matrix(runif(p),p,1)
#' y <- X %*% b + matrix(abs(rnorm(n)),n,1)
#' 
#' ## Setup mixture example
#' n <- 1e3
#' p <- 10
#' nnm <- generate_nnm(n,p)
#' set.seed(124)
#' X <- nnm$X
#' b <- double(p)
#' nComponents <- 3
#' k <- sample(1:p,nComponents,replace=FALSE)
#' b[k] <- matrix(runif(nComponents),ncol=1)
#' y <- X%*%b + 0.25*matrix(abs(rnorm(n)),n,1)
#' 
#' # Obtain solution to mixture problem
#' nnm_sol <- nnls_mm(y,X,runif(p))
#' 
nnls_mm <- function(y,X,b,max_iter=1e2,tol=1e-4) {
  # add checks for nonnegativity
  W <- apply(X,2,FUN=function(x) {return(x/sum(x))})
  b_last <- b
  for (iter in 1:max_iter) {
    b <- (t(W) %*% (y / (X%*%b_last))) * b_last
    if (norm(as.matrix(b - b_last),'f') < tol*(1 + norm(as.matrix(b_last), 'f')))
      break
    b_last <- b
  }
  return(list(b=b,iter=iter))
}

#' MM Algorithm - Plot NNM
#' 
#' \code{plot_nnm} Function for plotting nnm
#' 
#' @param nnm NNM object from generate_nnm function
#' 
#' @export
#' 
#' @examples
#' # Generate nonnegative matrix
#' n <- 1e3
#' p <- 10
#' nnm <- generate_nnm(n,p)
#' 
#' # Plot nonnegative matrix
#' plot_nnm(nnm)
#'   
#' @author Jocelyn T. Chi
#' 
plot_nnm <- function(nnm){
  x = values = ind = NULL
  a <- nnm$X
  nnm_stack <- stack(as.data.frame(a))
  nnm_stack$x <- rep(seq_len(nrow(a)), ncol(a))
  p <- qplot(x, values, data = nnm_stack, group = ind, colour = ind, geom = "line")
  p + theme_bw(base_size=14) + xlab("Frequency") + ylab("Intensity") + theme(legend.position = "none")
}

#' MM Algorithm - Plot NNM Objective
#' 
#' \code{plot_nnm_obj} Function for plotting the NNM Objective Function
#' 
#' @param y Nonnegative response
#' @param X Nonnegative design matrix
#' @param b Nonnegative initial regression vector
#' @param max_iter (Optional) Maximum number of iterations
#' 
#' @export
#' 
#' @examples
#' set.seed(12345)
#' n <- 100
#' p <- 3
#' X <- matrix(rexp(n*p,rate=1),n,p)
#' b <- matrix(runif(p),p,1)
#' y <- X %*% b + matrix(abs(rnorm(n)),n,1)
#' 
#' plot_nnm_obj(y,X,b)
#'   
#' @author Jocelyn T. Chi
#' 
plot_nnm_obj <- function(y,X,b,max_iter=100){
  bhat <- b
  loss <- double(max_iter)
  for (i in 1:max_iter) {
    bhat <- nnls_mm(y,X,bhat,max_iter=1)$b
    loss[i] <- 0.5*norm(as.matrix(y - X%*%bhat),'f')**2
  }
  x <- data.frame(1:max_iter)
  loss <- data.frame(loss)
  dat <- data.frame(x,loss)
  colnames(dat) <- c('x','loss')
  p <- ggplot(dat, aes(x=x,y=loss))
  p + geom_line() + theme_bw(base_size=14) + xlab("Iterates") + ylab("Value of the loss function")
}

#' MM Algorithm - Plotting the Spectroscopic Signal
#' 
#' \code{plot_spect} Function for plotting the spectroscopic signal
#' 
#' @param n Number of samples
#' @param nnm NNM object from generate_nnm function
#' @param y Nonnegative response
#' @param X Nonnegative design matrix
#' @param b Nonnegative initial regression vector
#' 
#' @export
#' 
#' @examples
#' # Setup mixture example
#' n <- 1e3
#' p <- 10
#' nnm <- generate_nnm(n,p)
#' 
#' set.seed(12345)
#' X <- nnm$X
#' b <- double(p)
#' nComponents <- 3
#' k <- sample(1:p,nComponents,replace=FALSE)
#' b[k] <- matrix(runif(nComponents),ncol=1)
#' y <- X%*%b + 0.25*matrix(abs(rnorm(n)),n,1)
#' 
#' plot_spect(n,y,X,b,nnm)
#'   
#' @author Jocelyn T. Chi
#' 
plot_spect <- function(n,y,X,b,nnm){
  t <- data.frame(nnm$t)
  y <- data.frame(y)
  dat <- data.frame(t,y)
  colnames(dat) <- c('t','y')
  p <- ggplot(dat, aes(x=t,y=y))
  p + geom_line() + theme_bw(base_size=14) + xlab("Frequency") + ylab("Intensity")
}

#' MM Algorithm - Plotting the Reconstruction
#' 
#' \code{plot_nnm_reconstruction} Function for plotting the nnm_sol reconstruction
#' 
#' @param nnm NNM object from generate_nnm function
#' @param X Nonnegative design matrix
#' @param nnm_sol Solution object from nnls_mm function
#' 
#' @export
#' 
#' @examples
#' # Setup mixture example
#' n <- 1e3
#' p <- 10
#' nnm <- generate_nnm(n,p)
#' 
#' set.seed(12345)
#' X <- nnm$X
#' b <- double(p)
#' nComponents <- 3
#' k <- sample(1:p,nComponents,replace=FALSE)
#' b[k] <- matrix(runif(nComponents),ncol=1)
#' y <- X%*%b + 0.25*matrix(abs(rnorm(n)),n,1)
#' 
#' # Obtain solution to mixture problem
#' nnm_sol <- nnls_mm(y,X,runif(p))
#' 
#' # Plot the reconstruction
#' plot_nnm_reconstruction(nnm,X,nnm_sol)
#' 
plot_nnm_reconstruction <- function(nnm,X,nnm_sol){
  x = NULL
  t <- data.frame(nnm$t)
  Xb <- data.frame(X%*%nnm_sol$b)
  dat <- data.frame(t,Xb)
  colnames(dat) <- c('t','Xb')
  p <- ggplot(dat, aes(x=t, y=Xb))
  p + geom_line() + theme_bw(base_size=14) + xlab("Frequency") + ylab("Reconstructed Intensity")
}

#' MM Algorithm - Plotting the True Signal
#' 
#' \code{plot_nnm_truth} Function for plotting the true mixture signal
#' 
#' @param nnm NNM object from generate_nnm function
#' @param X Nonnegative design matrix
#' @param b Nonnegative initial regression vector
#' 
#' @export
#' 
#' @examples
#' # Setup mixture example
#' n <- 1e3
#' p <- 10
#' nnm <- generate_nnm(n,p)
#' 
#' set.seed(12345)
#' X <- nnm$X
#' b <- double(p)
#' nComponents <- 3
#' k <- sample(1:p,nComponents,replace=FALSE)
#' b[k] <- matrix(runif(nComponents),ncol=1)
#' y <- X%*%b + 0.25*matrix(abs(rnorm(n)),n,1)
#' 
#' # Plot the truth
#' plot_nnm_truth(X,b,nnm)
#' 
plot_nnm_truth <- function(X,b,nnm){
  t <- data.frame(nnm$t)
  Xb <- data.frame(X%*%b)
  dat <- data.frame(t,Xb)
  colnames(dat) <- c('t','Xb')
  p <- ggplot(dat, aes(x=t, y=Xb))
  p + geom_line() + theme_bw(base_size=14) + xlab("Frequency") + ylab("True Intensity")
}

#' MM Algorithm - Plotting the NNMLS regression coefficients
#' 
#' \code{plot_nnm_coef} Function for plotting the NNMLS regression coefficients
#' 
#' @param nnm_sol Solution object from nnls_mm function
#' 
#' @export
#' 
#' @examples
#' # Setup mixture example
#' n <- 1e3
#' p <- 10
#' nnm <- generate_nnm(n,p)
#' 
#' set.seed(12345)
#' X <- nnm$X
#' b <- double(p)
#' nComponents <- 3
#' k <- sample(1:p,nComponents,replace=FALSE)
#' b[k] <- matrix(runif(nComponents),ncol=1)
#' y <- X%*%b + 0.25*matrix(abs(rnorm(n)),n,1)
#' 
#' # Obtain solution to mixture problem
#' nnm_sol <- nnls_mm(y,X,runif(p))
#' 
#' # Plot the regression coefficients
#' plot_nnm_coef(nnm_sol)
#' 
plot_nnm_coef <- function(nnm_sol){
  b <- data.frame(nnm_sol$b)
  x <- data.frame(1:nrow(b))
  dat <- data.frame(x,b)
  colnames(dat) <- c('x','b')
  p <- ggplot(dat, aes(x=x, y=b))
  p + geom_point() + theme_bw(base_size=14) + xlab("k") + ylab(expression(paste(b[k])))
}