
#---------------------------------------------------------------------
# cv wrapper for l2boost
# K-fold cross-validation for mean-squared prediction error
# makes use of mclapply
#---------------------------------------------------------------------
#' @title K-fold cross-validation using \code{\link{l2boost}}.
#'  
#' @description Calculate the K-fold cross-validation prediction error for \code{\link{l2boost}} models.
#' The prediction error is calculated using mean squared error (MSE). The optimal boosting step (\emph{m=opt.step})
#' is obtained by selecting the step \emph{m} resulting in the minimal MSE. 
#' 
#' @details
#' The cross-validation method splits the test data set into K mutually exclusive subsets. An \code{\link{l2boost}} model 
#' is built on K different training data sets, each created from a subsample of the full data set by sequentially leaving 
#' out one of the K subsets. The prediction error estimate is calculated by averaging the mean square error of each K test 
#' sets of the all of the K training datasets. The optimal step \emph{m} is obtained at the step with a minimal averaged 
#' mean square error.
#' 
#' The full \code{\link{l2boost}} model is run after the cross-validation models, on the full dateset. This model is
#' run for the full number of iteration steps \emph{M} and returned in the cv.l2boost$fit object. 
#' 
#' \code{\link{cv.l2boost}} only optimizes along the iteration count \emph{m} for a given value of \emph{nu}. This is
#' equivalent to an L1-regularization optimization. In order to optimize an elasticBoost model on the L2-regularization 
#' parameter lambda, a manual two way cross-validation can be obtained by sequentially optimizing over a range of lambda 
#' values, and selecting the lambda/opt.step pair resulting in the minimal cross-validated mean square error. See the 
#' examples below.
#' 
#' \code{\link{cv.l2boost}} uses the parallel package internally to speed up the cross-validation process on multicore 
#' machines. Parallel is packaged with base R >= 2.14, for earlier releases the multicore package provides the same 
#' functionality. By default, \code{\link{cv.l2boost}} will use all cores available except 1. Each fold is run on it's
#' own core and results are combined automatically. The number of cores can be overridden using the \emph{cores} function
#' argument. 
#' 
#' @param x the design matrix
#' @param y the response vector
#' @param K number of cross-validation folds (default: 10)
#' @param M the total number of iterations passed to \code{\link{l2boost}}.
#' @param nu l1 shrinkage paramater (default: 1e-4)
#' @param lambda l2 shrinkage parameter for elasticBoost (default: NULL = no l2-regularization)
#' @param type Type of l2boost fit with (default: discrete) see \code{\link{l2boost}} for description.
#' @param cores number of cores to parallel the cv analysis. If not specified, detects the 
#' number of cores. If more than 1 core, use n-1 for cross-validation. Implemented using 
#' multicore (mclapply), or clusterApply on Windows machines. 
#' @param trace Show computation/debugging output? (default: FALSE)
#' @param ... Additional arguments passed to \code{\link{l2boost}}
#' 
#' @seealso \code{\link{l2boost}}, \code{\link{plot.l2boost}}, 
#' \code{\link{predict.l2boost}} and \code{\link{coef.l2boost}}
#'
#' @return A list of cross-validation results:
#'  \item{call}{the matched call.}   
#'  \item{type}{Choice of l2boost algorithm from  "discrete", "hybrid", "friedman","lars". see \code{\link{l2boost}}}       
#'  \item{names}{design matrix column names used in the model}
#'  \item{nu}{The L1 boosting shrinkage parameter value}
#'  \item{lambda}{The L2 elasticBoost shrinkage parameter value}                   
#'  \item{K}{number of folds used for cross-validation}             
#'  \item{mse}{Optimal cross-validation mean square error estimate }      
#'  \item{mse.list}{list of \emph{K} vectors of mean square errors at each step \emph{m}}            
#'  \item{coef}{beta coefficient estimates from the full model at opt.step}     
#'  \item{coef.stand}{standardized beta coefficient estimates from full model at opt.step}    
#'  \item{opt.step}{optimal step m calculated by minimizing cross-validation error among all K training sets}           
#'  \item{opt.norm}{L1 norm of beta coefficients at opt.step}        
#'  \item{fit}{\code{\link{l2boost}} fit of full model}
#'  \item{yhat}{estimate of response from full model at opt.step}      
#' 
#' @examples
#' \dontrun{
#' #--------------------------------------------------------------------------
#' # Example: ElasticBoost simulation
#' # Compare l2boost and elasticNetBoosting using 10-fold CV
#' # 
#' # Elastic net simulation, see Zou H. and Hastie T. Regularization and 
#' # variable selection via the elastic net. J. Royal Statist. Soc. B, 
#' # 67(2):301-320, 2005
#' set.seed(1025)
#' dta <- elasticNetSim(n=100)
#' 
#' # The default values set up the signal on 3 groups of 5 variables,
#' # Color the signal variables red, others are grey.
#' sig <- c(rep("red", 15), rep("grey", 40-15))
#' 
#' # Set the boosting parameters
#' Mtarget = 1000
#' nuTarget = 1.e-2
#' 
#' # For CRAN, only use 2 cores in the CV method
#' cvCores=2
#' 
#' # 10 fold l2boost CV  
#' cv.obj <- cv.l2boost(dta$x,dta$y,M=Mtarget, nu=nuTarget, cores=cvCores)
#' 
#' # Plot the results
#' par(mfrow=c(2,3))
#' plot(cv.obj)
#' abline(v=cv.obj$opt.step, lty=2, col="grey")
#' plot(cv.obj$fit, type="coef", ylab=expression(beta[i]))
#' abline(v=cv.obj$opt.step, lty=2, col="grey")
#' plot(coef(cv.obj$fit, m=cv.obj$opt.step), cex=.5, 
#'   ylab=expression(beta[i]), xlab="Column Index", ylim=c(0,140), col=sig)
#'
#' # elasticBoost l1-regularization parameter lambda=0.1 
#' # 5 fold elasticNet CV
#' cv.eBoost <- cv.l2boost(dta$x,dta$y,M=Mtarget, K=5, nu=nuTarget, lambda=.1, cores=cvCores) 
#' 
#' # plot the results
#' plot(cv.eBoost)
#' abline(v=cv.eBoost$opt.step, lty=2, col="grey")
#' plot(cv.eBoost$fit, type="coef", ylab=expression(beta[i]))
#' abline(v=cv.eBoost$opt.step, lty=2, col="grey")
#' plot(coef(cv.eBoost$fit, m=cv.obj$opt.step), cex=.5, 
#'   ylab=expression(beta[i]), xlab="Column Index", ylim=c(0,140), col=sig)
#' }
#' @export cv.l2boost
#' @importFrom parallel mclapply
cv.l2boost <- function(x, y, K = 10, M = NULL, nu = 1e-4, lambda = NULL, trace = FALSE, 
                       type = c( "discrete", "hybrid", "friedman","lars"), cores=NULL,
                       ...) {
  call<-match.call()
  n<- length(y)
  
  # test for multicores
  if(is.null(cores)){
    num <- detectCores()
    # If the user doesn't know, then leave 1 core for OS type stuff. So the system
    # remains responsive.
    if(num > 1){
      cores = num-1
    }else{
      cores = 1
    }
  }else{
    if(cores > detectCores()){
      # it is ok to request more cores than the machine has, and can improve performance in
      # some instances. However, I'll still drop a warning.
      warning(paste("This CV analysis attempts to use more cores (", cores, ") than available (", 
                    detectCores(), "). \n !! This is not an error !!, analysis continuing."))
    } 
  }
  if(trace) cat("Number of cores:", cores, "\n")
  # Test for error conditions
  if(K > n) stop(paste("Number of folds (K=", K
                       , ") less than the number of observations in the design matrix (n=",
                       n, ")"))
  # define the folds
  # last fold is the full data and corresponds to the "primary object"
  all.folds <- split(sample(1:n), rep(1:K, length = n))
  all.folds[[K+1]] <- length(y) + 1
  
  if(cores ==1){
    # Single core function.
    eval.fold.obj <-lapply(1:length(all.folds), function(k) {
      eval.fold(k, K = K, all.folds = all.folds, x = x, y = y, M = M, nu = nu,
                lambda = lambda, trace = trace, type = type, ...=...)})
  }else{
    
    if (Sys.info()[1] == "Windows"){
      ## clusterApply() for Windows
      #       cl <- makeCluster(cores)
      #       eval.fold.obj <-clusterApply(cl=cl, k=1:length(all.folds), function(k) {
      #         eval.fold(k, K = K, all.folds = all.folds, x = x, y = y, M = M, nu = nu,
      #                   lambda = lambda, trace = trace, type = type, ...=...)})
      #       stopCluster(cl) # Don't forget to do this--I frequently do
      #      
      ## clusterApply() fails the build_win script, let's try it in serial.
      eval.fold.obj <-lapply(1:length(all.folds), function(k) {
        eval.fold(k, K = K, all.folds = all.folds, x = x, y = y, M = M, nu = nu,
                  lambda = lambda, trace = trace, type = type, ...=...)})
      
      # mclapply() for everybody else
    } else {
      # multicore cv-call
      eval.fold.obj <- mclapply(1:length(all.folds), function(k) {
        eval.fold(k, K = K, all.folds = all.folds, x = x, y = y, M = M, nu = nu,
                  lambda = lambda, trace = trace, type = type, ...=...)}, mc.cores=cores)
    }
  }
  #print(eval.fold.obj)
  # extract the mse
  mse.list <- lapply(1:K, function(k) {eval.fold.obj[[k]]$mse})
  M.step <- max(sapply(1:K, function(k){length(mse.list[[k]])}), na.rm = TRUE)
  # rewrite the mse list in a matrix format more conducive for parsing 
  cv.all <-  matrix(NA, nrow = M.step, ncol = K)
  for (k in 1:K) {
    cv.all[1:length(mse.list[[k]]), k] <- mse.list[[k]]
  }
  # determine the optimal number of steps
  mse.avg <- apply(cv.all, 1, mean, na.rm = TRUE)
  opt.step <- min(c(which(mse.avg == min(mse.avg, na.rm = TRUE)), M))
  mse <- mse.avg[opt.step]
  # prediction using full data for the optimal number of steps
  fit.all <- eval.fold.obj[[K+1]]$obj
  yhat <- predict.l2boost(fit.all)$yhat.path[[opt.step]]
  pred <- predict.l2boost(fit.all, type = "coef")
  opt.norm <- sum(abs(pred$coef.path[[opt.step]]), na.rm = TRUE)
  # add names to the optimal coefficient path
  opt.coef.path <- pred$coef.path[[opt.step]]
  opt.coef.stand.path <- pred$coef.stand.path[[opt.step]]
  names(opt.coef.path) <- names(opt.coef.stand.path) <- fit.all$names
  
  # return the object
  object <- list(call=call, type=type,nu=nu,K=K,lambda=lambda, 
                 fit = fit.all,
                 mse = mse,
                 mse.list = mse.list,
                 yhat = yhat,
                 coef = opt.coef.path,
                 coef.stand = opt.coef.stand.path,
                 opt.step = opt.step,
                 opt.norm = opt.norm,
                 names = fit.all$names)
  class(object) <- c("l2boost", "cv")
  invisible(object)
}
