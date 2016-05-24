#' Elastic Multinomial Logistic Regression
#'
#' This function identifies a multinomial logistic regression model with
#' phase-variablity using elastic methods
#'
#' @param f matrix (\eqn{N} x \eqn{M}) of \eqn{M} functions with \eqn{N} samples
#' @param y vector of size \eqn{M} labels {1,2,...,m} for m classes
#' @param time vector of size \eqn{N} describing the sample points
#' @param B matrix defining basis functions (default = NULL)
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
#' \item{Loss}{logistic loss}
#' \item{type}{model type ('mlogistic')}
#' @keywords srvf alignment regression
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Elastic Functional Logistic Regression with Application to Physiological Signal Classification,
#'  Electronic Journal of Statistics (2014), submitted.
elastic.mlogistic <- function(f, y, time, B=NULL, df=20, max_itr=20,
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

  # Code Labels
  m = max(y)
  Y = matrix(0,N,m)
  for (ii in 1:N){
    Y[ii, y[ii]] = 1
  }

  if (smooth_data){
    f = smooth.data(f,sparam)
  }

  # Create B-Spline Basis if none provided
  if (is.null(B)){
    B = bs(time, df=df, degree=4, intercept=TRUE)
  }
  Nb = ncol(B)

  # Compute q-function of the functional data
  tmp = gradient.spline(f,binsize,smooth_data)
  f = tmp$f
  q = tmp$g/sqrt(abs(tmp$g)+eps)

  gam = kronecker(matrix(1,1,N),seq(0,1,length.out=M))

  itr = 1
  LL = rep(0, max_itr)
  while(itr <= max_itr) {
    cat(sprintf("Iteration: r=%d\n", itr))
    # align data
    fn = matrix(0, M, N)
    qn = matrix(0, M, N)
    for (ii in 1:N){
      fn[,ii] = approx(time,f[,ii],xout=(time[length(time)]-time[1])*gam[,ii] + time[1])$y
      qn[,ii] = f_to_srvf(fn[,ii], time)
    }

    # Find alpha and beta using l_bfgs
    Phi = matrix(1, N, Nb+1)
    for (ii in 1:N){
      for (jj in 2:(Nb+1)){
        Phi[ii,jj] = trapz(time,qn[,ii] * B[, jj-1])
      }
    }
    b0 = rep(0,m*(Nb+1))
    out = optim(b0, mlogit_loss, gr = mlogit_gradient, Phi, Y,
                method = "L-BFGS-B", control = list(maxit=200,pgtol=1e-10))
    b = out$par

    B0 = array(b,c(Nb+1, m))
    alpha = B0[1,]
    beta = matrix(0,M,m)
    for (ii in 1:m){
      beta[,ii] = B %*% B0[2:(Nb+1),ii]
    }

    # compute the Loss
    LL[itr] = mlogit_loss(b,Phi,Y)

    # find gamma
    k=1
    gamma_new<-foreach(k = 1:N, .combine=cbind,.packages="fdasrvf") %dopar% {
      gam = mlogit_warp(alpha, beta, time, q[,k], Y[k,])
    }

    if (pvecnorm(gam-gamma_new,2) < 1e-5){
      break
    }else{
      gam = gamma_new
    }

    itr = itr + 1
  }
  gam = gamma_new

  if (parallel){
    stopCluster(cl)
  }

  return(list(alpha=alpha, beta=beta, fn=fn, qn=qn, gamma=gam, q=q, B=B,
              b=b[2:length(b)], Loss=LL[1:itr], n_classes=m, mode='mlogistic')  )
}
