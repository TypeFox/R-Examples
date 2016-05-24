#' @title Thresholded EEBoost
#' 
#' @description 
#' Run the thresholded EEBoost procedure.
#' 
#' @export 
#' 
#' @param Y Vector of outcomes.
#' @param X Matrix of predictors. Will be automatically scaled using the \code{\link{scale}} function.
#' @param EE.fn Estimating function taking arguments \code{Y}, \code{X}, and parameter vector \code{b}.
#' @param b.init Initial parameter values. For variable selection, typically start with a vector of zeroes (the default).
#' @param eps Step length. Default is 0.01, value should be relatively small.
#' @param maxit Maximum number of iterations. Default is 1000.
#' @param itertrack Indicates whether or not diagnostic information should be printed out at each iteration. Default is \code{FALSE}.
#' @param reportinterval If \code{itertrack} is \code{TRUE}, how many iterations the algorithm should wait between each diagnostic report.
#' @param stop.rule Rule for stopping the iterations before \code{maxit} is reached. Possible values are \code{"on.repeat"} and \code{"pct.change"}.
#' See 'Details' for more information.
#' @param thresh Threshold parameter for ThrEEBoost.
#' 
#' @return A matrix with \code{maxit} rows and \code{ncol(X)} columns, with each row containing the parameter vector from an iteration of ThrEEBoost.
#' 
#' @details 
#' \code{threeboost} Implements a thresholded version of the EEBoost algorithm described in \emph{Wolfson (2011, JASA)}.
#' EEBoost is a general-purpose method for variable selection which can be applied whenever inference would be based on an estimating equation.
#' The package currently implements variable selection based on the Generalized Estimating Equations, but can also accommodate 
#' user-provided estimating functions. Thresholded EEBoost is a generalization which allows multiple variables to enter the model at each boosting step. 
#' Thresholded EEBoost with thresholding parameter = 1 is equivalent to EEBoost.
#' 
#' Typically, the boosting procedure is run for \code{maxit} iterations, producing \code{maxit} models defined by a set of regression coefficients.
#' An additional step (e.g. model scoring, cross-validated estimate of prediction error) is needed to select a final model. However, an alternative is to stop the iterations
#' before \code{maxit} is reached. The user can request this feature by setting \code{stop.rule} to one of the following options:
#' 
#' \itemize{
#' \item \code{"on.repeat"}: Sometimes, ThrEEBoost will alternate between stepping on the same two directions, usually indicating numerical problems. Setting \code{stop.rule="on.oscillate"} will terminate the algorithm if this happens.
#' \item \code{"pct.change"}: Stop if, for conseuctive iterations, the sum of the magnitudes of the elements of the estimating equation changes by < 1\%. 
#' }
#' 
#' @examples
#' library(Matrix)
#' 
#' # Generate some test data - uses 'mvtnorm' package
#' n <- 30
#' n.var <- 50
#' clust.size <- 4
#' B <- c(rep(2,5),rep(0.2,5),rep(0.05,10),rep(0,n.var-20))
#' mn.X <- rep(0,n.var)
#' sd.X <- 0.5
#' rho.X <- 0.3
#' cov.sig.X <- sd.X^2*((1-rho.X)*diag(rep(1,10)) + rho.X*matrix(data=1,nrow=10,ncol=10))
#' sig.X <- as.matrix( Matrix::bdiag(lapply(1:(n.var/10),function(x) { cov.sig.X } ) ) )
#' sd.Y <- 0.5
#' rho.Y <- 0.3
#' indiv.Sig <- sd.Y^2*( (1-rho.Y)*diag(rep(1,4)) + rho.Y*matrix(data=1,nrow=4,ncol=4) )
#' sig.list <- list(length=n)
#' for(i in 1:n) { sig.list[[i]] <- indiv.Sig }
#' Sig <- Matrix::bdiag(sig.list)
#' indiv.index <- rep(1:n,each=clust.size)
#' sig.Y <- as.matrix(Sig)
#' if(require(mvtnorm)) {
#' X <- mvtnorm::rmvnorm(n*clust.size,mean=mn.X,sigma=sig.X)
#' mn.Y <- X %*% B
#' ## Correlated continuous outcome
#' Y <- mvtnorm::rmvnorm(1,mean=mn.Y,sigma=sig.Y) 
#' } else { stop('Need mvtnorm package to generate correlated example data.') }
#' 
#' ## Define the Gaussian GEE estimating function with independence working correlation
#' mu.Lin <- function(eta){eta}
#' g.Lin <- function(m){m}
#' v.Lin <- function(eta){rep(1,length(eta))}
#'
#'  EE.fn.ind <- function(Y,X,b) { 
#'  ee.GEE(Y,X,b,
#'  mu.Y=mu.Lin,
#'  g.Y=g.Lin,
#'  v.Y=v.Lin,
#'  aux=function(...) { ee.GEE.aux(...,mu.Y=mu.Lin,g.Y=g.Lin,v.Y=v.Lin) },
#'  id=indiv.index,
#'  corstr="ind")
#' }
#' 
#' ## These two give the same result
#' coef.mat <- eeboost(Y,X,EE.fn.ind,maxit=250)
#' coef.mat2 <- geeboost(Y,X,id=indiv.index,family="gaussian",corstr="ind",maxit=250)$coefmat
#' 
#' par(mfrow=c(1,2))
#' coef_traceplot(coef.mat)
#' coef_traceplot(coef.mat2)
#' 
#' @seealso \code{\link{geeboost}} for an example of how to call (Thr)EEBoost with a custom estimating function.
#' 
#' Wolfson, J. \href{http://pubs.amstat.org/doi/abs/10.1198/jasa.2011.tm10098}{EEBoost: A general method for prediction and variable selection using estimating equations.} Journal of the American Statistical Association, 2011.
#' 
threeboost <- function(Y,X,EE.fn,b.init=rep(0,ncol(X)),eps=0.01,maxit=1000,itertrack=FALSE,reportinterval=1,stop.rule="on.repeat",thresh=1) {
  
  Y <- as.vector(Y)
  X <- as.matrix(X)
  
  if(thresh < 0 | thresh > 1) { stop('ERROR: Threshold must be between 0 and 1.')}
  
  b.new <- b.old <- b.init
  ee.val <- rep(0,length(b.new))
  it <- 1
  
  B <- matrix(data=NA,nrow=maxit,ncol=length(b.init))
  scale.X <- scale(X)
  
  if(any(is.nan(scale.X))) { 
    stop("Covariate matrix cannot be scaled due to degenerate covariate.")
  }
  
  while(it <= maxit) {
    b.old2 <- b.old
    b.old <- b.new
    ee.val.old <- ee.val
    
    ee.val <- round(EE.fn(Y,scale.X,b.old),8) ## To counter numerical problems in the next step
    three.val <- (abs(ee.val)>=thresh*max(abs(ee.val))) ## Threshold the estimating equation
    b.new <- b.old + eps*sign(ee.val)*three.val
    
    if(itertrack==TRUE & (it %% reportinterval==0)) {
      cat(paste("----------------\n Iteration ",it,"\n Stepping on direction ",which.max(abs(ee.val)),"\n Max abs element of EE vector: ",max(abs(ee.val)),"\n"))
      print(rbind(b.new[which(b.new!=0)],which(b.new!=0)))
    }
    if(stop.rule=="on.repeat" & ( all(b.new==b.old2) | all(b.new==b.old))) {
      print(paste("Stopped due to oscillating values at iteration",it))
      B[it,] <- b.new
      return(B[1:it,])
    }
    if(stop.rule=="pct.change" & (sum(abs(ee.val))-sum(abs(ee.val.old)))/sum(abs(ee.val.old)) < 0.01) {
      print(paste("Stopped due to < 1% change in sum of absolute EE elements at iteration",it))
      B[it,] <- b.new
      return(B[1:it,])
    }
    B[it,] <- b.new
    it <- it+1
  }
  B
}