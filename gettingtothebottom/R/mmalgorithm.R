#' MM Algorithm - Softhreshold Function
#' 
#' \code{softhreshold} Function for computing the softhreshold
#' 
#' @param x Vector of values to be softhresholded
#' @param lambda Softhreshold parameter
#' 
#' @export
#' 
#' @examples
#' x <- seq(-10,10,1)
#' softhreshold(x,lambda=3)
#'   
#' @author Jocelyn T. Chi
#' 
softhreshold <- function(x,lambda) {
  x_gelambda <- which(x >= lambda)
  x_lelambda <- which(x <= -lambda)
  sol <- double(length(x))
  sol[x_gelambda] <- x[x_gelambda] - lambda
  sol[x_lelambda] <- x[x_lelambda] + lambda
  return(sol)
}

#' MM Algorithm - Make Z
#' 
#' \code{makeZ} Function for making the Z matrix
#' 
#' @param M Matrix containing observed entries
#' @param lambda Softhreshold parameter
#' 
#' @export
#' 
#' @examples
#' A <- matrix(rnorm(9),3,3)
#' makeZ(A,lambda=3)
#'   
#' @author Jocelyn T. Chi
#' 
makeZ <- function(M,lambda){
  m <- dim(M)[1]
  n <- dim(M)[2]
  M.svd <- svd(M)
  r <- min(m,n)
  U <- M.svd$u[,1:r,drop=FALSE]
  D <- M.svd$d[1:r]
  V <- M.svd$v[,1:r,drop=FALSE]
  S <- diag(softhreshold(D,lambda),r,r)
  Z <- U %*% S %*% t(V)
return(Z)
}

#' MM Algorithm - Make Y
#' 
#' \code{makeY} Function for making the Y matrix
#' 
#' @param X Matrix containing observed entries
#' @param Z Matrix containing last iterates
#' @param omega Vector containing indices of unobserved entries (by column)
#' 
#' @import Matrix
#' 
#' @export
#' 
#' @examples
#' n <- 5
#' A <- matrix(rnorm(n^2),n,n)
#' omega <- c(1,5,8,10,16,23)
#' Z <- Matrix(0,n,n,sparse=TRUE)
#' makeY(A,Z,omega)
#'   
#' @author Jocelyn T. Chi
#' 
makeY <- function(X,Z,omega){
  m <- dim(X)[1]
  n <- dim(X)[2]
  Y <- matrix(0,m,n)
  omega_compl <- setdiff(1:(m*n),omega)
  Y[omega_compl] <- X[omega_compl]
  Y[omega] <- Z[omega]
  return(Y)
}

#' MM Algorithm - Matrix Completion
#' 
#' \code{matrixcomplete} Function for performing matrix completion using a majorization-minimization algorithm given data matrix X
#' 
#' @param X Data matrix to be completed
#' @param Z Matrix containing last iterates
#' @param omega Vector containing indices of unobserved entries
#' @param lambda Softhreshold parameter
#' @param maxiter (Optional) Max number of iterations (Default: 100)
#' @param tol (Optional) Tolerance for convergence (Default: 1e-4)
#' @param liveupdates (Optional) If FALSE, no notification will be given upon completion of each iteration. (Default: TRUE)
#' 
#' @export
#' 
#' @examples
#' # (Examples not run)
#' # Generate an m-by-n test matrix of rank r
#' # seed <- 12345
#' # m <- 1000
#' # n <- 1000
#' # r <- 5
#' # T <- testmatrix(m,n,r,seed=seed)
#' 
#' # Add some noise to the test matrix
#' # E <- 0.1*matrix(rnorm(m*n),m,n)
#' # A <- T + E
#' 
#' # Obtain a vector of unobserved entries
#' # temp <- makeOmega(m,n,percent=0.5)
#' # omega <- temp$omega
#' 
#' # Remove unobserved entries from test matrix
#' # X <- A
#' # X[omega] <- NA
#' 
#' # Make initial model matrix Z and find initial lambda
#' # Z <- matrix(0,m,n)
#' # lambda <- init.lambda(X,omega)
#' 
#' # Example (Not run)
#' # Sol <- matrixcomplete(X,Z,omega,lambda)
#'   
#' @author Jocelyn T. Chi
#'
matrixcomplete <- function(X,Z,omega,lambda,maxiter=100,tol=1e-4,liveupdates=TRUE){
  loss <- double(maxiter)
  for (i in 1:maxiter){
    Zold <- as.matrix(Z)
    Y <- makeY(X,Zold,omega)
    Z <- makeZ(Y,lambda) 
    # compute relative change in Z
    difference <- diff_norm(Z,Zold,omega)
    loss[[i]] <- (1/2)*norm(Y-Z,'F')^2 + lambda*sum(svd(Z)$d)
    if (difference < tol*(1 + norm(as.matrix(Zold),'f'))) {
      cat('Optimization completed.')
      cat('\n')
      break
    }
    else{
      if (liveupdates == TRUE){
        cat(c('Completed iteration ',i,'.'))
        cat('\n')
      }
    }
  }
  return(list(sol=Z,iter=i,loss=loss[1:i]))
}

#' MM Algorithm - Generate Test Matrix
#' 
#' \code{testmatrix} Function for generating random rank-r matrix
#' 
#' @param r Rank of matrix to be generated (r >= 2)
#' @param m Number of rows in matrix to be generated
#' @param n Number of columns in matrix to be generated
#' @param seed Random seed
#' 
#' @export
#' 
#' @examples
#' m <- 100
#' n <- 1000
#' r <- 5
#' testmatrix(m,n,r)
#'  
#' @author Jocelyn T. Chi
#'
testmatrix <- function(m,n,r,seed=123){
#  x <- c(1:5)
  set.seed(seed)
  M <- matrix(rnorm(n*m),m,n)
  M.svd <- svd(M)
  U <- M.svd$u
  D <- M.svd$d
  V <- M.svd$v
  Mat <- U[,1:r] %*% diag(D[1:r]) %*% t(V[,1:r])
  return(Mat)
}

#' MM Algorithm - Generate Omega
#' 
#' \code{makeOmega} Function for generating omega vector of missing values
#' 
#' @param m Number of rows in matrix to be generated
#' @param n Number of columns in matrix to be generated
#' @param percent Percent missing in matrix
#' @param seed Random seed
#' 
#' @export
#' 
#' @examples
#' m <- 1000
#' n <- 1000
#' percent <- 0.75
#' omega <- makeOmega(m,n,percent)
#'  
#' @author Jocelyn T. Chi
#'
makeOmega <- function(m,n,percent,seed=123){
  if (percent > 1 | percent < 0){
    cat('The percentage of missing values must range between 0 and 1.')
    break
  }
  set.seed(seed)
  tot <- m*n
  vec <- 1:tot
  Omega <- sample(vec,round(percent*tot),replace=FALSE) 
  OmegaComplement <- setdiff(vec,Omega)
  return(list(omega=Omega,omegac=OmegaComplement))
}

#' MM Algorithm - Plot the Softhreshold Function
#' 
#' \code{plot_softhreshold} Function for plotting the softhreshold function
#' 
#' @param from The starting value of the sequence of inputs into the function
#' @param to The ending value of the sequence of inputs into the function
#' @param lambda Lambda value for softhreshold function
#' 
#' @import ggplot2
#' 
#' @export
#' 
#' @examples
#' plot_softhreshold(-5,5,3)
#'  
#' @author Jocelyn T. Chi
#'
plot_softhreshold <- function(from=-5, to=5,lambda=3) {
  if (lambda < from | lambda > to) {
    print('Please select a lambda between the starting and end values of the inputs')
    break
  }
  x <- seq(from,to,(to-from)/100)
  f <- softhreshold(x,lambda=lambda)
  x <- data.frame(x)
  f <- data.frame(f)
  dat <- data.frame(x,f)
  p <- ggplot(dat, aes(x=x,y=f)) + ylim(from,to) + xlim(from,to)
  p + geom_line(colour="blue") + geom_abline(colour = "gray") + theme_bw(base_size=14) + xlab("x") + ylab("Softhreshold value")
}

#' MM Algorithm - Normed Difference
#' 
#' \code{diff_norm} Function for finding the normed difference between two matrices based on vector containing indices of differing elements
#' 
#' @param X Original data matrix
#' @param Z Model matrix for comparison
#' @param omega Set of unobserved indices
#' 
#' @export
#' 
#' @examples
#' Z <- matrix(rnorm(9,0,1),3,3)
#' X <- matrix(rnorm(9,0,2),3,3)
#' omega <- c(2,5,6)
#' 
#' diff_norm(X,Z,omega)
#'  
#' @author Jocelyn T. Chi
#'
diff_norm <- function(X,Z,omega){
  x <- X[omega]
  z <- Z[omega]
  diffn <- norm(as.matrix(x-z),'F')  
  return(diffn)
}

#' MM Algorithm - Find the best fit lambda for a given problem based on an initial guess for lambda
#' 
#' \code{solutionpaths} Function for finding the best fit lambda for a given problem based on an initial guess for lambda
#' 
#' @param A Original data matrix (no unobserved entries)
#' @param X Data matrix (with unobserved entries)
#' @param Z Initial model matrix
#' @param omega Vector of unobserved entries in the data matrix X
#' @param lambda.start Initial value for lambda
#' @param tol (Optional) Tolerance for convergence (Default: 1e-4)
#' @param liveupdates (Optional) Set to TRUE to view progress of comparisons. (Default: FALSE)
#' @param lambdaseq_length (Optional) Length of lambda sequence for convergence. (Default: 20)
#' 
#' @export
#' 
#' @examples
#' # Generate a test matrix
#' seed <- 12345
#' m <- 100
#' n <- 100
#' r <- 3
#' T <- testmatrix(m,n,r,seed=seed)
#' 
#' # Add some noise to the test matrix
#' E <- 0.1*matrix(rnorm(m*n),m,n)
#' A <- T + E
#' 
#' # Obtain a vector of unobserved entries
#' temp <- makeOmega(m,n,percent=0.5)
#' omega <- temp$omega
#' 
#' # Remove unobserved entries from test matrix
#' X <- A
#' X[omega] <- NA
#' 
#' # Make initial model matrix Z and find initial lambda
#' Z <- matrix(0,m,n)
#' lambda.start <- init.lambda(X,omega)
#' lambdaseq_length=20
#' tol <- 1e-2
#' 
#' ans <- solutionpaths(A,X,Z,omega,lambda.start,tol=tol,
#'    liveupdates=TRUE,lambdaseq_length=lambdaseq_length)
#' 
#' @author Jocelyn T. Chi
#'
solutionpaths <- function(A,X,Z,omega,lambda.start,tol=1e-4,liveupdates=FALSE,lambdaseq_length=20){
  L <- log10(lambda.start)
  lambdaseq <- makeLambdaseq(L,lambdaseq_length)
  results.lambda <- double(lambdaseq_length)
  results.difference <- double(lambdaseq_length)
  Zlist <- vector("list",lambdaseq_length)
  for (i in 1:lambdaseq_length) {    
    ans <- matrixcomplete(X,Z,omega,lambda=lambdaseq[i],maxiter=1e3,tol=tol,liveupdates=liveupdates)
    Sol <- ans$sol
    results.difference[i] <- diff_norm(Sol,A,omega)
    Zlist[[i]] <- Sol
    if (liveupdates==TRUE) {
      cat('Completed results for lambda ',i,' of ',length(lambdaseq),'.')
      cat('\n') 
    }
    Z <- Sol
  }
  return(list(results=Zlist,lambda=lambdaseq,error=results.difference))
}

#' MM Algorithm - Function for making sequence of lambdas for solution paths
#' 
#' \code{makeLambdaseq} Function for making sequence of lambdas for solution paths given starting lambda value and desired sequence length
#' 
#' @param L Initial lambda value
#' @param lambdaseq_length Desired length of lambda sequence
#' 
#' @export
#' 
#' @examples
#' makeLambdaseq(11,20)
#'  
#' @author Jocelyn T. Chi
#'
makeLambdaseq <- function(L,lambdaseq_length){

  lambdaseq <- 10^seq(from=L,to=L-2,length.out=lambdaseq_length)
  return(lambdaseq)
}

#' MM Algorithm - Function for plotting the imputed values against the truth for minimum error solution
#' 
#' \code{plot_solpaths_error} Function for plotting the imputed values against the truth for minimium error solution found using solutionpaths function
#' 
#' @param A Initial test matrix of fully observed entries
#' @param omega Vector of unobserved entries
#' @param ans Results from solpaths function
#' 
#' @export
#' 
#' @examples
#' # Generate a test matrix
#' seed <- 12345
#' m <- 100
#' n <- 100
#' r <- 3
#' T <- testmatrix(m,n,r,seed=seed)
#' 
#' # Add some noise to the test matrix
#' E <- 0.1*matrix(rnorm(m*n),m,n)
#' A <- T + E
#' 
#' # Obtain a vector of unobserved entries
#' temp <- makeOmega(m,n,percent=0.5)
#' omega <- temp$omega
#' 
#' # Remove unobserved entries from test matrix
#' X <- A
#' X[omega] <- NA
#' 
#' # Make initial model matrix Z and find initial lambda
#' Z <- matrix(0,m,n)
#' lambda.start <- init.lambda(X,omega)
#' lambdaseq_length=20
#' tol <- 1e-2
#' 
#' ans <- solutionpaths(A,X,Z,omega,lambda.start,tol=tol,
#'    liveupdates=FALSE,lambdaseq_length=lambdaseq_length)
#' 
#' plot_solpaths_error(A,omega,ans)
#'  
#' @author Jocelyn T. Chi
#'
plot_solpaths_error <- function(A, omega, ans){
  Truth = Imputation = NULL
  min_error <- which(ans$error == min(ans$error))
  dat <- data.frame(A[omega],ans$results[[min_error]][omega])
  colnames(dat) <- c('Truth','Imputation')
  p <- ggplot(dat, aes(x=Truth,y=Imputation))
  p + geom_point() + geom_abline(colour = "gray") + theme_bw(base_size=14) + xlab("True values") + ylab("Imputed values")
}

#' MM Algorithm - Plot results of solutionpaths function
#' 
#' \code{plot_solutionpaths} Function for plotting results of the solutionpaths function
#' 
#' @param results Results from the solutionpaths function
#' 
#' @export
#' 
#' @examples
#' # Generate a test matrix
#' seed <- 12345
#' m <- 100
#' n <- 100
#' r <- 3
#' T <- testmatrix(m,n,r,seed=seed)
#' 
#' # Add some noise to the test matrix
#' E <- 0.1*matrix(rnorm(m*n),m,n)
#' A <- T + E
#' 
#' # Obtain a vector of unobserved entries
#' temp <- makeOmega(m,n,percent=0.5)
#' omega <- temp$omega
#' 
#' # Remove unobserved entries from test matrix
#' X <- A
#' X[omega] <- NA
#' 
#' # Make initial model matrix Z and find initial lambda
#' Z <- matrix(0,m,n)
#' lambda.start <- init.lambda(X,omega)
#' lambdaseq_length=20
#' tol <- 1e-2
#' 
#' ans <- solutionpaths(A,X,Z,omega,lambda.start,tol=tol,
#'    liveupdates=FALSE,lambdaseq_length=lambdaseq_length)
#' 
#' # Plot using results from solutionpaths function
#' plot_solutionpaths(ans)
#'  
#' @author Jocelyn T. Chi
#'
plot_solutionpaths <- function(results) {
  logLambda = Error = rdLambda = NULL
  n <- length(results$lambda)
  lambdas <- data.frame(results$lambda)
  rounded_lambdas <- round(lambdas,2)
  error <- data.frame(results$error)
  loglambdas <- log10(lambdas)
  dat <- data.frame(cbind(loglambdas,lambdas,rounded_lambdas,error))
  colnames(dat) <- c('logLambda','Lambda','rdLambda','Error')
  min_error <- which(dat$Error == min(dat$Error))
  min_error_lambda <- dat$Lambda[min_error]
  vline <- log10(min_error_lambda)
  p <- ggplot(dat, aes(x=logLambda,y=Error))
  p <- p + geom_point() + geom_vline(xintercept = vline, colour="blue") + geom_text(aes(label=rdLambda, size=1), hjust=1, vjust=2, angle=30) + ylim(min(dat$Error)-5,max(dat$Error)+5)
  p + theme_bw(base_size=14) + theme(legend.position = "none") + xlab(expression(paste(log[10](lambda),"(rounded values beneath points)"))) + ylab("Error")
}

#' MM Algorithm - Initial lambda
#' 
#' \code{init.lambda} Function for finding an initial value for lambda
#' 
#' @param X Original data matrix
#' @param omega Set of unobserved indices
#' 
#' @export
#' 
#' @examples
#' # Generate a test matrix
#' seed <- 12345
#' m <- 100
#' n <- 100
#' r <- 3
#' T <- testmatrix(m,n,r,seed=seed)
#' 
#' # Add some noise to the test matrix
#' E <- 0.1*matrix(rnorm(m*n),m,n)
#' A <- T + E
#' 
#' # Obtain a vector of unobserved entries
#' temp <- makeOmega(m,n,percent=0.5)
#' omega <- temp$omega
#' 
#' # Remove unobserved entries from test matrix
#' X <- A
#' X[omega] <- NA
#' 
#' init.lambda(X,omega)
#'  
#' @author Jocelyn T. Chi
#'
init.lambda <- function(X,omega){
  X[omega] <- 0
  X.svd <- svd(X)
  lambda <- X.svd$d[1]
  return(lambda)
}
