#' Williams' Type Test Statistic.
#' 
#' @description
#' Calculates a Williams' type test statistic for a constrained linear mixed effects model.
#' 
#' @rdname w.stat
#' 
#' 
#' @param theta estimated coefficients.
#' @param cov.theta covariance matrix of the (unconstrained) coefficients. 
#' @param B matrix to obtain the global contrast. 
#' @param A   matrix of linear constraints.
#' @param ... additional arguments, to enable custom test statistic functions.
#' 
#' @details
#' See \code{\link{create.constraints}} for an example of \code{A}. Argument \code{B} is similar, but defines the global contrast for a Williams' type test statistic. This is the largest hypothesized difference in the constrained coefficients. So for an increasing simple order, the test statistic is the difference between the two extreme coefficients, \eqn{\theta_1}{theta_1} and \eqn{\theta_{p_1}}{theta_p1}, divided by the standard error (unconstrained). For an umbrella order order, two contrasts are considered, \eqn{\theta_1}{theta_1} to \eqn{\theta_{s}}{theta_s}, and \eqn{\theta_{p_1}}{theta_p1} to \eqn{\theta_{s}}{theta_s}, each divided by the appropriate unconstrained standard error. A general way to express this statistic is:
#' 
#' \deqn{W = max \theta_{B[i,2]} - \theta_{B[i,1]} / sqrt( VAR( \theta_{B[i,2]} - \theta_{B[i,1]} ) )}{W = max theta_{B[i,2]} - theta_{B[i,1]} / sqrt( VAR( theta_{B[i,2]} - theta_{B[i,1]} ) )}
#' 
#' where the numerator is the difference in the constrained estimates, and the standard error in the denominator is based on the covariance matrix of the unconstrained estimates.
#' 
#' The function \code{w.stat.ind} does the same, but uses the \code{A} matrix which defines all of the individual constraints, and returns a test statistic for each constraints instead of taking the maximum.
#'
#' @return
#' Output is a numeric value.
#' 
#' @note
#' See \code{\link{lrt.stat}} for information on creating custom test statistics.
#' 
#' @examples
#' theta  <- exp(1:4/4)
#' th.cov <- diag(4)
#' X1     <- matrix( 0 , nrow=1 , ncol=4 )
#' const  <- create.constraints( P1=4 , constraints=list(order='simple' ,
#'                                                     decreasing=FALSE) )
#' 
#' w.stat( theta , th.cov , const$B , const$A )
#' 
#' w.stat.ind( theta , th.cov , const$B , const$A )
#' 
#' @export
#' 

##
## Williams' type statistic (global)
##
w.stat <- function( theta , cov.theta , B , A , ...  ){
  
  stats <- vector( "numeric",length=nrow(B) )
  ctd   <- diag( cov.theta )
  
  stats <- apply( B , 1 , 
                  FUN=function(b,theta,cov,ctd){
                    std <- sqrt( ctd[b[1]] + ctd[b[2]] - 2*cov.theta[b[1],b[2]] )
                    (theta[b[2]]-theta[b[1]])/std
                  }, theta=theta, cov=cov.theta, ctd=ctd)
  
  test.stat <- max( stats )
  tnames    <- names(theta)
  ii        <- min( which( stats == max(stats) ) )  
  
  if( !is.null(tnames) ){
    names(test.stat) <- paste0( tnames[B[ii,2]] , " - ", tnames[B[ii,1]] )
  } else{
    names(test.stat) <- paste0( "Theta ", B[ii,2] , " - Theta ", B[ii,1] )
  }
  
  return( test.stat )
  
}


#' Williams' type statistic (individual)
#' 
#' @rdname w.stat
#' @export
#' 

w.stat.ind <- function( theta , cov.theta , B , A , ...  ){
  
  stats <- vector( "numeric",length=nrow(A) )
  ctd   <- diag( cov.theta )
  
  stats <- apply( A , 1 , 
                  FUN=function(a,theta,cov,ctd){
                    std <- sqrt( ctd[a[1]] + ctd[a[2]] - 2*cov.theta[a[1],a[2]] )
                    (theta[a[2]]-theta[a[1]])/std
                  }, theta=theta, cov=cov.theta, ctd=ctd)
  
  return( stats )
  
}



#' Likelihood ratio type statistic (global)
#' 
#' @description
#' Calculates the likeihood ratio type test statistic (under Normality assumption) for a 
#' constrained linear mixed effects model. This is the default test statistic for \pkg{CLME}.
#' 
#' 
#' @param theta estimated coefficients.
#' @param theta.null coefficients estimated under the null hypothesis.
#' @param cov.theta covariance matrix of the (unconstrained) coefficients. 
#' @param ... additional arguments, to enable custom test statistic functions.
#' 
#' @return
#' Output is a numeric value.
#' 
#' @note
#' This is an internal function, unlikely to be useful outside of \link{CLME-package}. To define custom functions, the arguments available are:
#' 
#' \code{theta}, \code{theta.null}, \code{cov.theta}, \code{B}, \code{A}, \code{Y}, \code{X1}, \code{X2}, \code{U}, \code{tsq}, \code{ssq}, \code{Nks}, and \code{Qs}.
#' 
#' Of the additional arguments, \code{B} and \code{A} are identical to those produced by \code{\link{create.constraints}}. The rest, \code{Y}, \code{X1}, \code{X2}, \code{U}, \code{tsq}, , \code{ssq}, \code{Nks}, and \code{Qs}, are equivalent to arguments to \code{\link{clme_em}}.
#' 
#' Custom functions must produce numeric output. Output may have length greater than 1, which corresponds to testing multiple global hypotheses.
#' 
#' @seealso
#' \code{\link{clme_em}},
#' \code{\link{w.stat}}
#' 
#' 
#' @examples
#' data( rat.blood )
#' cons <- list(order = "simple", decreasing = FALSE, node = 1 )
#' 
#' clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
#'                  constraints = cons, seed = 42, nsim = 0)
#' 
#' # Individually compute lrt statistic
#' lrt.stat(clme.out$theta, clme.out$theta.null, clme.out$cov.theta )
#' 
#' @export
#' 
#' 



lrt.stat <- function( theta, theta.null, cov.theta, ... ){  
  theta.diff <- theta - theta.null  
  test.stat  <- c( t(theta.diff) %*% cov.theta %*% theta.diff )
  
  names(test.stat) <- "Bootstrap LRT"
  
  # Return test statistic
  return(test.stat)
  
}


