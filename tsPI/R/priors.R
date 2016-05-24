#' Compute different types of importance weights based on Jeffreys's prior
#'
#' These functions compute different types of importance weights based on Jeffreys's priors used in \code{\link{arima_pi}}.
#'
#' @export
#' @rdname priors
#' @name jeffreys
#' @useDynLib tsPI, .registration=TRUE
#' @seealso \code{\link{arima_pi}}.
#' @param psi vector containing the ar and ma parameters (in that order).
#' @param xreg matrix or data frame containing the exogenous variables
#' (not including the intercept which is always included for non-differenced series)
#' @param p number of ar parameters
#' @param q number of ma parameters
#' @param n length of the time series
approx_joint_jeffreys<-function(psi, xreg = NULL, p, q, n){

  # Large sample approximation for Jeffreys joint prior p(beta,sigma,psi)=p(psi)/sigma=
  # sqrt(|XinvVX'||I(psi)|), where I is the large sample approximation of information matrix
  # Note that I_sigma,psi->0 when n->infinity so that term is omitted
  # (equals to the case where they are assumed independent a priori)
  #
  phi <- if (p > 0) psi[1:p] else NULL
  theta <- if (q > 0) psi[p + 1:q] else NULL
  chol_im <- try(chol(information_arma(phi, theta)), silent = TRUE)
  if (class(chol_im) != "try-error") {
    sqrt_det <- prod(diag(chol_im))
    if (!is.null(xreg)){
      sqrt_det * sqrt(det(t(cbind(1, xreg)) %*% solve(toeplitz(acv_arma(phi, theta, n = n)),cbind(1,xreg))))
    } else sqrt_det*sqrt(sum(solve(toeplitz(acv_arma(phi, theta, n = n)), matrix(1,nrow = n))))
  }  else 0
}
#' @export
#' @rdname priors
#'
approx_marginal_jeffreys<-function(psi, p, q){
  # Large sample approximation for Jeffreys marginal prior p(psi) =
  # sqrt(|I(psi)|), where I is the large sample approximation of information matrix
  phi <- if (p > 0) psi[1:p] else NULL
  theta <- if (q > 0) psi[p + 1:q] else NULL
  chol_im <- try(chol(information_arma(phi, theta)), silent = TRUE)
  if (class(chol_im) != "try-error") {
    prod(diag(chol_im))
  } else 0
}
#' @export
#' @rdname priors
exact_joint_jeffreys<-function(psi, xreg = NULL, p, q, n){
  # Exact Jeffreys's joint prior
  phi <- if(p>0) psi[1:p] else NULL
  theta <- if(q>0) psi[(p+1):(p+q)] else NULL
  V<-toeplitz(acv_arma(phi, theta,n = n))
  invV<-solve(V)
  dV<-dacv_arma(phi, theta,n = n)
  informationMatrixPsi<-matrix(0,p+q,p+q)
  informationMatrixPsiSigma<-numeric(p+q)
  dVMatrices<-vector("list",length=p+q)
  dx<-1 + 0:(n - 1) * (n + 1)
  for(i in 1:(p+q)){
    dVMatrices[[i]]<-toeplitz(dV[,i])
    informationMatrixPsiSigma[i]<-sum((invV %*% dVMatrices[[i]])[dx])
  }
  for(i in 1:(p+q))
    for(j in 1:(p+q))
      informationMatrixPsi[i,j]<-sum((invV %*% dVMatrices[[i]]%*%invV %*% dVMatrices[[j]])[dx])

  if(!is.null(xreg)){
    sqrt(det(t(cbind(1,xreg))%*%invV%*%cbind(1,xreg))*
        det(informationMatrixPsi-1/(2*n)*informationMatrixPsiSigma%*%informationMatrixPsiSigma))
  } else sqrt(matrix(1,ncol=n)%*%invV%*%matrix(1,nrow=n)*
      det(informationMatrixPsi-1/(2*n)*informationMatrixPsiSigma%o%informationMatrixPsiSigma))
}
#' @export
#' @rdname priors
exact_marginal_jeffreys <- function(psi, p , q, n){
  # Exact Jeffreys's marginal prior
  phi <- if(p>0) psi[1:p] else NULL
  theta <- if(q>0) psi[(p+1):(p+q)] else NULL
  V<-toeplitz(acv_arma(phi, theta, n = n))
  invV<-solve(V)
  dV<-dacv_arma(phi, theta,n = n)
  informationMatrixPsi<-matrix(0,p+q,p+q)
  dVMatrices<-vector("list",length=p+q)
  dx<-1 + 0:(n - 1) * (n + 1)
  for(i in 1:(p+q))
    dVMatrices[[i]]<-toeplitz(dV[,i])
  for(i in 1:(p+q))
    for(j in 1:(p+q))
      informationMatrixPsi[i,j]<-sum((invV %*% dVMatrices[[i]]%*%invV %*% dVMatrices[[j]])[dx])
  sqrt(det(informationMatrixPsi))
}

