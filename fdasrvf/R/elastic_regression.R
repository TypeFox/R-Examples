#' Elastic Linear Regression
#'
#' This function identifies a regression model with phase-variablity
#' using elastic methods
#'
#' @param f matrix (\eqn{N} x \eqn{M}) of \eqn{M} functions with \eqn{N} samples
#' @param y vector of size \eqn{M} responses
#' @param time vector of size \eqn{N} describing the sample points
#' @param B matrix defining basis functions (default = NULL)
#' @param lam scalar regularization parameter (default=0)
#' @param df scalar controling degrees of freedom if B=NULL (default=20)
#' @param max_itr scalar number of iterations (default=20)
#' @param smooth_data smooth data using box filter (default = F)
#' @param sparam number of times to apply box filter (default = 25)
#' @param parallel enable parallel mode using \code{\link{foreach}} and
#'   \code{doParallel} pacakge
#' @param cores set number of cores to use with \code{doParallel} (default = 2)
#' @return Returns a list containing
#' \item{alpha}{model intercept}
#' \item{beta}{regressor function}
#' \item{fn}{aligned functions - matrix (\eqn{N} x \eqn{M}) of \eqn{M} functions with \eqn{N} samples}
#' \item{qn}{aligned srvfs - similar structure to fn}
#' \item{gamma}{warping functions - similar structure to fn}
#' \item{q}{original srvf - similar structure to fn}
#' \item{B}{basis matrix}
#' \item{b}{basis coefficients}
#' \item{SSE}{sum of squared errors}
#' \item{type}{model type ('linear')}
#' @keywords srvf alignment regression
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Elastic Functional Logistic Regression with Application to Physiological Signal Classification,
#'  Electronic Journal of Statistics (2014), submitted.
elastic.regression <- function(f, y, time, B=NULL, lam=0, df=20, max_itr=20,
                               smooth_data = FALSE, sparam = 25, parallel = FALSE,
                               cores=2){

  if (parallel){
    cl = makeCluster(cores)
    registerDoParallel(cl)
  } else
  {
    registerDoSEQ()
  }

  binsize = mean(diff(time))
  eps = .Machine$double.eps
  M = nrow(f)
  N = ncol(f)
  f0 = f

  if (smooth_data){
    f = smooth.data(f,sparam)
  }

  # Create B-Spline Basis if none provided
  if (is.null(B)){
    B = bs(time, df=df, degree=4, intercept=TRUE)
  }
  Nb = ncol(B)

  # second derivative
  Bdiff = matrix(0,M,Nb)
  for (ii in 1:Nb){
    Bdiff[,ii] = gradient(gradient(B[,ii],binsize), binsize)
  }

  # Compute q-function of the functional data
  tmp = gradient.spline(f,binsize,smooth_data)
  f = tmp$f
  q = tmp$g/sqrt(abs(tmp$g)+eps)

  gam = kronecker(matrix(1,1,N),seq(0,1,length.out=M))

  itr = 1
  SSE = rep(0, max_itr)
  while(itr <= max_itr) {
    cat(sprintf("Iteration: r=%d\n", itr))
    # align data
    fn = matrix(0, M, N)
    qn = matrix(0, M, N)
    for (ii in 1:N){
      fn[,ii] = approx(time,f[,ii],xout=(time[length(time)]-time[1])*gam[,ii] + time[1])$y
      qn[,ii] = f_to_srvf(fn[,ii], time)
    }

    # OLS using basis
    Phi = matrix(1, N, Nb+1)
    for (ii in 1:N){
      for (jj in 2:(Nb+1)){
        Phi[ii,jj] = trapz(time,qn[,ii] * B[, jj-1])
      }
    }

    R = matrix(0, Nb+1, Nb+1)
    for (ii in 2:Nb+1){
      for (jj in 2:(Nb+1)){
        R[ii,jj] = trapz(time,Bdiff[,ii-1] * Bdiff[,jj-1])
      }
    }

    xx = t(Phi) %*% Phi
    inv_xx = solve(xx + lam * R)
    xy = t(Phi) %*% y
    b = inv_xx %*% xy

    alpha = b[1]
    beta = B %*% b[2:(Nb+1)]

    # compute the SSE
    int_X = rep(0, N)
    for (ii in 1:N){
      int_X[ii] = trapz(time, qn[,ii] * beta)
    }

    SSE[itr] = sum((y-alpha-int_X)^2)

    # find gamma
    gamma_new<-foreach(k = 1:N, .combine=cbind,.packages="fdasrvf") %dopar% {
      gam = regression_warp(beta, time, q[,k], y[k], alpha)
    }

    if (pvecnorm(gam-gamma_new,2) < 1e-5){
      break
    }else{
      gam = gamma_new
    }

    itr = itr + 1
  }

  # last step with centering of gam
  gamI = SqrtMeanInverse(t(gam))
  gamI_dev = gradient(gamI, 1/(M-1))
  beta = approx(time,beta,xout=(time[length(time)]-time[1])*gamI +
                      time[1])$y*sqrt(gamI_dev)

  for (k in 1:N){
    qn[,k] = approx(time,qn[,k],xout=(time[length(time)]-time[1])*gamI +
                         time[1])$y*sqrt(gamI_dev)
    fn[,k] = approx(time,fn[,k],xout=(time[length(time)]-time[1])*gamI +
                         time[1])$y
    gam[,k] = approx(time,gam[,k],xout=(time[length(time)]-time[1])*gamI +
                       time[1])$y
  }

  if (parallel){
    stopCluster(cl)
  }

  return(list(alpha=alpha, beta=beta, fn=fn, qn=qn, gamma=gam, q=q, B=B,
              b=b[2:length(b)], SSE=SSE[1:itr], mode='linear')  )
}
