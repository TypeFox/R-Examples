#' @title Construct MInt object
#' 
#' @description Constructs a MInt object that maintains the data and parameter estimates for
#' the underlying Poisson-multivariate normal hierarchical model.
#' 
#' @param y A file path to the response matrix.
#' @param x A file path to the design matrix.
#' @param fmla An object of class ``\code{\link{formula}}'' (or one that can be
#'   coerced to that class): a symbolic description of the model to be fitted.
#'   
#' @return mint A MInt object.
#'   
#' 
#' @examples
#' x <- system.file("extdata", "x.txt", package="MInt");
#' y <- system.file("extdata", "y.txt", package="MInt");
#' m <- mint(y,x,fmla = ~feature1 + feature2)
#' 
#' @import testthat
#' @export
mint <- function(y, x, fmla= ~1) {
  # Constructor
  
  # Implement basic parameter checks.
  
  # Construct the object.
  mfit <- structure(list(data=list(design=x, response=y, fmla=fmla), 
                         param=list(), optim=list()), class="mint");
  
  mfit$data$ready <- FALSE;
  mfit$optim$ready <- FALSE;
  mfit$param$ready <- FALSE;
  
  return(mfit);
}




## Function implementations

#' @title Estimate parameters
#' 
#' @description
#' This function performs iterative conditional modes to obtain
#' maximum \emph{a posteriori} estimates for \eqn{\beta} (covariate coefficients), 
#' \eqn{w} (latent abundances), and \eqn{P} (the precision matrix).
#' 
#' 
#' @param mfit - a MInt model object.
#' 
#' @return A MInt model object with the following attributes:
#'  \item{optim}{List containing optimization details}
#'  \item{optim$lambda}{Value of the L1 penalty used during optimization}
#'  \item{data}{List containing the raw data}
#'  \item{data$design}{File path of the design matrix}
#'  \item{data$response}{File path of the response matrix}
#'  \item{data$fmla}{Formula used to model each response in terms of the design variables}
#'  \item{data$y}{Raw numerical data for the response matrix}
#'  \item{data$xd}{Design matrix in categorical form}
#'  \item{data$x}{Design matrix in numerical form}
#'  \item{param}{List containing parameter estimates}
#'  \item{param$beta}{p-covariates x o-responses matrix of regression coefficients}
#'  \item{param$w}{n-samples x o-responses matrix of latent abundances}
#'  \item{param$P}{o-responses x o-responses precision matrix}
#'  
#' @examples
#' x <- system.file("extdata", "x.txt", package="MInt");
#' y <- system.file("extdata", "y.txt", package="MInt");
#' m <- mint(y,x,fmla = ~feature1 + feature2)
#' m <- estimate(m)
#'  
#'  
#' @import glasso
#' @import trust
#' @import MASS
#' @export
estimate <- function(mfit) {
  
  
  # Read data, initialize optimizer, and initialize parameters if 
  # not already done.
  if (!mfit$data$ready) mfit <- read_data(mfit);
  if (!mfit$optim$ready) mfit <- initialize_optimizer(mfit);
  if (!mfit$param$ready) mfit <- initialize_parameters(mfit);
  
  while (!mfit$optim$converged & mfit$optim$itr <= 1000) {
    mfit <- update_beta(mfit);
    mfit <- update_w(mfit);
    mfit <- update_mu(mfit);
    mfit <- update_P(mfit);
    
    mfit <- check_convergence(mfit);
    
    cat("Iteration: ", mfit$optim$itr, 
        " max(deltaP): ", sprintf("%0.6f", max(mfit$optim$deltaP)), 
        " Objective: ", mfit$optim$objective_trace[length(mfit$optim$objective_trace)], 
        " Converged: ", mfit$optim$converged, "\n");
  }
  
  return(mfit);
}

#' @title Bootstrap
#'   
#' @description This function bootstraps a model learned by \code{estimate} to
#' obtain confidence intervals on each parameter.
#' 
#' @param mfit A MInt model object.
#' @param nboot The number of bootstraps to perform.
#' @param seed Random number generator seed.
#'   
#' @return A MInt object.
#'   
#'   We should export this at some point.
bootstrap <- function(mfit,nboot=10,seed=1) {
  
  set.seed(seed);
  bs <- vector("list", nboot)
  for (i in 1 : nboot){
    mfit_copy <- mfit;
    
    # Resample
    bsidx <- sample.int(nrow(mfit$data$y), nrow(mfit$data$y), TRUE)
    mfit_copy$data$y <- mfit$data$y[bsidx,];
    mfit_copy$data$x <- as.matrix(mfit$data$x[bsidx,]);
    mfit_copy$data$xd <- as.data.frame(mfit$data$xd[bsidx,], row.names=rownames(mfit$data$xd), col.names=colnames(mfit$data$xd));
    
    # Reinitialize
    mfit_copy <- initialize_optimizer(mfit_copy);
    mfit_copy <- initialize_parameters(mfit_copy);
    
    # Estimateq
    mfit_copy <- estimate(mfit_copy);
    
    bs[[i]] <- mfit_copy$param;
    
    cat("Bootstrap replicate ", i, " complete.\n");
  }
  mfit$param$bootstrap <- bs;
  return(mfit);
}


read_data <- function(mfit) {
  # Read in the data.
  # xd  maintains the categorical representation of the data.
  # x   maintains a complete numerical representation of the data. 
  #     Continuous variables are kept as is, whereas categorical variables are
  #     split into their binary representations.
  # y   is the response matrix. 
  
  mfit <- read_response(mfit);
  mfit <- read_design(mfit);
  
  mfit$data$ready <- TRUE;
  return(mfit)
}

initialize_optimizer <- function(mfit) {
  if (is.null(mfit$data$y)) {
    stop("Uninitialized response matrix.")
  }
  
  o <- ncol(mfit$data$y);
  n <- nrow(mfit$data$y);
  
  mfit$optim$objective_constant <- -sum(lfactorial(mfit$data$y)) - o*n*log(2*pi)/2 + (o^2)*log(n/4);
  mfit$optim$tolerance <- 1e-4;
  mfit$optim$lambda <- 0.01;
  mfit$optim$itr <- 0;
  
  # Parameter clamps (for debugging)
  mfit$optim$clamp$P <- FALSE;
  mfit$optim$clamp$mu <- FALSE;
  mfit$optim$clamp$w <- FALSE;
  mfit$optim$clamp$beta <- FALSE;
  
  mfit$optim$converged <- FALSE;
  mfit$optim$ready <- TRUE;
  
  return(mfit);
}

initialize_parameters <- function(mfit) {
  # Initialize the parameters.
  # beta
  # Start with a log-normal fit. 
  mfit$param$xb <- matrix(, nrow=nrow(mfit$data$y), ncol=ncol(mfit$data$y));
  mfit <- initialize_beta(mfit);
  mfit <- xb(mfit);
  
  # w
  # Set to be the residuals of the log-normal fit.
  mfit$param$w <- log(mfit$data$y+1) - mfit$param$xb;
  
  # mu
  mfit <- update_mu(mfit);
  
  # S and P
  # Set to be cov(w) and pseudoinverse of cov(w)
  mfit$param$S <- pop_cov(mfit$param$w, mfit$param$mu);
  mfit$param$P <- MASS::ginv(mfit$param$S); 
  rownames(mfit$param$P) <- colnames(mfit$data$y);
  colnames(mfit$param$P) <- colnames(mfit$data$y);
  rownames(mfit$param$S) <- colnames(mfit$data$y);
  colnames(mfit$param$S) <- colnames(mfit$data$y);
  mfit$optim$oldP <- mfit$param$P; # Store a copy of P for convergence checking.
  
  # Update hyperparameter lambda to its conditional mode.
  mfit$optim$lambda <- 2*(length(mfit$param$mu)^2)/( nrow(mfit$data$y)*sum(abs(mfit$param$P)) );
  
  mfit$param$ready <- TRUE;
  return(mfit)
}




## Update functions
update_beta <- function(mfit) {
  
  if (!mfit$optim$clamp$beta) {
    for (j in 1 : ncol(mfit$data$y)) {
      y <- mc(mfit$data$y,j); # extract a response vector
      w <- mc(mfit$param$w,j); # extract the corresponding latent abundance vector
      fmla <- update(mfit$data$fmla, y ~ . + offset(w));
      
      # Assemble the data and run the GLM.
      data <- cbind(y, mfit$data$xd[[j]], w);
      fit <- glm(fmla, family=poisson(), data=data);
      
      mfit$param$beta[,j] <- fit$coefficients;
    }
    
    mfit <- xb(mfit);
  }
  
  return(mfit);
}

update_mu <- function(mfit) {
  if (!mfit$optim$clamp$mu){
    mfit$param$mu <- matrix(replicate(ncol(mfit$data$y),0),nrow=1); #matrix(colMeans(mfit$param$w),nrow=1);
    colnames(mfit$param$mu) <- colnames(mfit$data$y);
  }
  
  return(mfit)
}

update_P <- function(mfit) {
  if (!mfit$optim$clamp$P){
    res <- glasso::glasso(pop_cov(mfit$param$w, mfit$param$mu), mfit$optim$lambda);
    #start="warm", w.init=mfit$param$S, wi.init=mfit$param$P);
    mfit$param$P <- res$wi;
    mfit$param$S <- res$w;
    
    rownames(mfit$param$P) <- colnames(mfit$data$y);
    colnames(mfit$param$P) <- colnames(mfit$data$y);
  }
  
  return(mfit);
}

update_w <- function(mfit) {
  if (!mfit$optim$clamp$w){
    for (ridx in 1:nrow(mfit$param$w)){
      wopt <- trust::trust(objfun=wfun, parinit=mfit$param$w[ridx,], 
                    rinit=5, rmax=100, minimize=FALSE, ridx=ridx, y=mfit$data$y, param=mfit$param);
      mfit$param$w[ridx,] <- wopt$argument; # the maximizer.
    }
  }
  
  return(mfit);
}

check_convergence <- function(mfit) {
  mfit <- evaluate_objective(mfit);
  mfit$optim$itr <- mfit$optim$itr + 1;
  mfit$optim$deltaP <- abs(mfit$param$P - mfit$optim$oldP); 
  
  mfit$optim$converged <- max(mfit$optim$deltaP) < mfit$optim$tolerance 
  
  mfit$optim$oldP <- mfit$param$P;
  
  return(mfit);
}



## Accessory functions
initialize_beta <- function(mfit) {
  # Compute this as a loop in case there are response-specific 
  # design matrices (i.e. mfit$data$x[[i]] != mfit$data$x[[j]], for i != j)
  mfit$param$beta <- matrix(, nrow=ncol(mfit$data$x[[1]]), ncol=ncol(mfit$data$y));
  for (j in 1 : ncol(mfit$data$y)) {
    mfit$param$beta[,j] <- MASS::ginv(mfit$data$x[[j]]) %*% log(mfit$data$y[,j] + 1);
  }
  rownames(mfit$param$beta) <- colnames(mfit$data$x[[1]]);
  colnames(mfit$param$beta) <- colnames(mfit$data$y);
  return(mfit);
}

xb <- function(mfit) {
  # Assumes mfit$param$xb is already preallocated.
  # Compute this as a loop in case there are response-specific 
  # design matrices (i.e. mfit$data$x[[i]] != mfit$data$x[[j]], for i != j)
  for (j in 1 : ncol(mfit$data$y)) {
    mfit$param$xb[,j] <- mfit$data$x[[j]] %*% mc(mfit$param$beta,j);
  }
  dimnames(mfit$param$xb) <- dimnames(mfit$data$y);
  return(mfit);
}

read_response <- function(mfit, row.names="Observations") {
  mfit$data$y <- as.matrix(read.table(mfit$data$response, header=TRUE, row.names=row.names));
  return(mfit);
}

read_design <- function(mfit, row.names="Observations") {
  # Check if the supplied 'design' input is a file or a directory.
  # If it's a directory, then there should be separate files of response 
  # specific design matrices.
  #   - Design matrix files for each response should be labeled as
  #     <response_name>.txt, where <response_name> is an element of 
  #     colnames(mfit$data$y)
  # If it's a file, then this it should be a single design matrix that 
  # all response variables share.
  
  if (is.null(mfit$data$y)) {
    stop("Uninitialized response matrix. Unknown response variable names.");
  }
  
  mfit$data$xd <- list();
  mfit$data$x <- list();
  
  if (file.info(mfit$data$design)$isdir) {
    
    # Change directory to the design directory, load all of the design matrices in the order of
    # colnames(mfit$data$y), and change back to the original directory.
    old <- getwd();
    setwd(mfit$data$design)
    
    for (rn in colnames(mfit$data$y)){
      mfit$data$xd[[rn]] <- model.frame(mfit$data$fmla, read.table(paste0(rn, ".txt"), header=TRUE, row.names=row.names));
      mfit$data$x[[rn]] <- model.matrix(mfit$data$fmla, mfit$data$xd[[rn]]);
    }
    
    setwd(old);
    
  } else if (file.exists(mfit$data$design)) {
    
    # Load the design matrix and copy it for all responses in colnames(mfit$data$y)
    # We are copying it so that all functions deal with the same data structure regardless
    # of whether there is a single design matrix or one for each response.
    xd <- model.frame(mfit$data$fmla, read.table(mfit$data$design, header=TRUE, row.names=row.names));
    x <- model.matrix(mfit$data$fmla, xd);
    
    for (rn in colnames(mfit$data$y)){
      mfit$data$xd[[rn]] <- xd;
      mfit$data$x[[rn]] <- x;
    }
    
  } else {
    stop("Unrecognized path or file for design matrix.")
  }
  
  return(mfit);
}

last_char <- function(s) return(substr(s, nchar(s), nchar(s)))

is.defined <- function(x) !is.null(x)

mr <- function(m,idx) {
  # Single row slice of matrix m.
  t(as.matrix(m[idx,]));
}

mc <- function(m,idx) {
  # Single column slice of matrix m.
  as.matrix(m[,idx]);
}

pop_cov <- function(x, mu) {
  n <- nrow(x);
  xc <- scale(x,center=mu,scale=FALSE);
  S <- t(xc)%*%xc/n;
  return(S);
}

wfun <- function(wrow, ...){
  args <- list(...);
  
  # Compute the function value, gradient, and Hessian.
  eta <- exp( mr(args$param$xb, args$ridx) + wrow);
  muc <- wrow - args$param$mu;
  lf <- muc %*% args$param$P;
  
  # Return function values consistent with MAXIMIZATION.
  f <- sum( mr(args$y, args$ridx)*wrow - eta )  - lf%*%t(muc)/2;
  g <- t(mr(args$y, args$ridx) - eta - lf);
  H <- -args$param$P - diag(as.vector(eta));
  
  # Return 
  list(value=f, gradient=g, hessian=H);
}

evaluate_objective <- function(mfit) {
  n <- nrow(mfit$data$y);
  eta <- exp(mfit$param$xb + mfit$param$w);
  dv <- determinant(mfit$param$P,logarithm=TRUE);
  ldP <- dv$modulus;
  
  poisson_component <- sum(mfit$data$y*log(eta) - eta);
  gaussian_component <- n*ldP/2 - n*sum(pop_cov(mfit$param$w,mfit$param$mu)*mfit$param$P)/2;
  laplacian_component <- -n*mfit$optim$lambda*sum(abs(mfit$param$P))/2;
  
  if (is.null(mfit$optim$objective_trace)) mfit$optim$objective_trace <- vector(mode="numeric");
  mfit$optim$objective_trace[length(mfit$optim$objective_trace)+1] <- 
    poisson_component + gaussian_component + laplacian_component + mfit$optim$objective_constant;
  
  return(mfit);
}

