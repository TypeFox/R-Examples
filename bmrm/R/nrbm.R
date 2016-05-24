
#' Convex and non-convex risk minimization with L2 regularization and limited memory
#' 
#' Use algorithm of Do and Artieres, JMLR 2012 to find w minimizing: 
#' f(w) = 0.5*LAMBDA*l2norm(w) + riskFun(w)
#' where riskFun is either a convex or a non-convex risk function.
#' @param riskFun the risk function to use in the optimization (e.g.: hingeLoss, softMarginVectorLoss). 
#'   The function must evaluate the loss value and its gradient for a given point vector (w).
#' @param LAMBDA control the regularization strength in the optimization process. 
#'   This is the value used as coefficient of the regularization term.
#' @param MAX_ITER the maximum number of iteration to perform. 
#'   The function stop with a warning message if the number of iteration exceed this value
#' @param EPSILON_TOL control optimization stoping criteria: 
#'   the optimization end when the optimization gap is below this threshold
#' @param w0 initial weight vector where optimization start
#' @param maxCP mximal number of cutting plane to use to limit memory footprint
#' @param convexRisk a length 1 logical telling if the risk function riskFun is convex. 
#'    If TRUE, use CRBM algorithm; if FALSE use NRBM algorithm from Do and Artieres, JMLR 2012
#' @return the optimal weight vector (w)
#' @export
#' @import LowRankQP
#' @examples
#'   set.seed(123)
#'   X <- matrix(rnorm(4000*200), 4000, 200)
#'   beta <- c(rep(1,ncol(X)-4),0,0,0,0)
#'   Y <- X%*%beta + rnorm(nrow(X))
#'   w <- nrbm(ladRegressionLoss(X/100,Y/100),maxCP=50)
#'   layout(1)
#'   barplot(w)
nrbm <- function(riskFun,LAMBDA=1,MAX_ITER=1000L,EPSILON_TOL=0.01,w0=0,maxCP=100L,convexRisk=TRUE) {
  R <- riskFun(w0)
  at <- as.vector(gradient(R))
  w0 <- rep(w0,length.out=length(at))  
  bt <- as.vector(R) - crossprod(w0,at)
  
  A <- matrix(numeric(0),0L,length(at))
  b <- numeric(0)
  inactivity.score <- numeric(0)
  a0 <- b0 <- NULL
  s0 <- 0
  s <- numeric(0)
    
  # find minimizer of the underestimator function
  # and update aggregated cutting plane
  optimize <- function() {
    # add aggregated cutting cutting plane to A and b
    A <- rbind(a0,A)
    b <- c(b0,b)
    
    Ale <- matrix(1,1L,nrow(A)+1L)
    H <- matrix(0,1L+nrow(A),1L+nrow(A))
    H[-1,-1] <- tcrossprod(A)    
    opt <- LowRankQP(H,c(0,-LAMBDA*b),Ale,1,rep(1,nrow(A)+1L),method="LU")
    alpha <- opt$alpha[-1L]
    
    # update aggregated cutting plane
    inactivity.score <<- inactivity.score + pmax(1-alpha[-1L],0)
    a0 <<- colSums(alpha * A)
    b0 <<- sum(alpha * b)
    
    w <- as.vector(-crossprod(A,alpha) / LAMBDA)
    R <- max(0,A %*% w + b)
    return(list(w = w, obj = LAMBDA*0.5*crossprod(w)+R))
  }  

  opt <- list(w=w0)
  ub.w <- w0
  ub.R <- R
  ub <- LAMBDA*0.5*crossprod(w0) + R
  st <- 0
  is.newbest <- TRUE
  for (i in 1:MAX_ITER) {    
    # add the new cutting plane to the working set
    cp <- head(order(inactivity.score,na.last=FALSE),n=maxCP)
    A <- rbind(at,A[cp,])
    b <- c(bt,b[cp])
    s <- c(st,s[cp])
    if (is.newbest) {
      inactivity.score[is.na(inactivity.score)] <- 0
      inactivity.score<- c(NA_real_,inactivity.score[cp])
    } else {
      inactivity.score <- c(0,inactivity.score[cp])
    }
      
    # optimize the underestimator
    opt <- optimize()
    lb <- opt$obj
    
    # test for the end of convergence
    cat(sprintf("%d:%.3f %.3f\n",i,ub-lb,ub))
    if (ub-lb < EPSILON_TOL) break
    
    # estimate loss at the new underestimator optimum
    R <- riskFun(opt$w)
    f <- LAMBDA*0.5*crossprod(opt$w) + R
    
    # deduce parameters of the new cutting plane
    at <- as.vector(gradient(R))
    bt <- R - crossprod(opt$w,at)
    
    if (!convexRisk) {
      # solve possible conflicts with the new cutting plane
      if (f<ub) {
        st <- 0
        s <- s + 0.5*LAMBDA*crossprod(ub.w-opt$w)
        s0 <- s0 + 0.5*LAMBDA*crossprod(ub.w-opt$w)
        b <- pmin(b,R - (A %*% opt$w) - s)
        b0 <- pmin(b0,R - crossprod(a0,opt$w) - s0)
      } else { # null step
        st <- 0.5*LAMBDA*crossprod(opt$w-ub.w)
        if (ub.R < st + crossprod(at,ub.w) + bt) {
          U <- ub.R - crossprod(at,ub.w) - st
          L <- ub - crossprod(at,opt$w) - 0.5*LAMBDA*crossprod(opt$w)
          if (L<=U) {
            bt <- L
          } else {
            at <- -LAMBDA*ub.w
            bt <- ub - 0.5*LAMBDA*crossprod(opt$w) - crossprod(at,opt$w)
          }
        }
      }
    }
    is.newbest <- (f<ub)
    if (is.newbest) {
      ub <- f
      ub.w <- opt$w
      ub.R <- R
    }
  }
  if (i >= MAX_ITER) warning('max # of itertion exceeded')
  return(ub.w)
}







