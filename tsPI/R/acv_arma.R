#' Compute a theoretical autocovariance function of ARMA process
#'
#' Function \code{acv_arma} computes a theoretical autocovariance function of ARMA process.
#'
#' @export
#' @rdname acv_arma
#' @name acv_arma
#' @seealso \code{\link{dacv_arma}}.
#' @param phi vector containing the AR parameters
#' @param theta vector containing the MA parameters
#' @param n length of the time series
#' @return vector of length n containing the autocovariances
#' @examples
#'
#' ## Example from Brockwell & Davis (1991, page 92-94)
#' ## also in help page of ARMAacf (from stats)
#' n <- 0:9
#' answer <- 2^(-n) * (32/3 + 8 * n) /(32/3)
#' acv <- acv_arma(c(1.0, -0.25), 1.0, 10)
#' all.equal(acv/acv[1], answer)
#'
acv_arma<-  function(phi,theta,n){
  acv <- numeric(n)
  model <- SSMarima(phi,theta,n=1)
  Tt<-model$T
  prodtp<-model$P1
  acv[1]<- prodtp[1]
  for(i in 2:n){
    prodtp<-Tt%*%prodtp
    acv[i]<-prodtp[1]
  }
  acv
}

#' Compute the partial derivatives of theoretical autocovariance function of ARMA process
#'
#' Function \code{dacv_arma} computes the partial derivatives of theoretical autocovariance function of ARMA process
#'
#' @export
#' @rdname dacv_arma
#' @seealso \code{\link{acv_arma}}.
#
#' @param phi vector containing the AR parameters
#' @param theta vector containing the MA parameters
#' @param n length of the time series
#' @return matrix containing the partial derivatives autocovariances,
#' each column corresponding to one parameter of vector (phi,theta) (in that order)
dacv_arma<-function(phi,theta,n){
  model<-SSMarima(phi,theta,n=1)
  p<-length(phi)
  q<-length(theta)
  dV<-matrix(0,n,p+q)

  Tt<-model$T
  m<-dim(Tt)[1]
  m2<-m^2
  kroneckerTtTt<-diag(m2)-kronecker(Tt,Tt)
  if(p>0){
    dTt<-matrix(0,m,m)
    kroneckerTtdTt<-matrix(0,m2,m2)
    for(i in 1:p){
      dTt[]<-0
      kroneckerTtdTt[]<-0
      dTt[i,1]<-1
      for(j in 1:m)
        for(k in 1:m)
          kroneckerTtdTt[(j-1)*m+i,(k-1)*m+1]<-Tt[j,k]

      kroneckerTtdTt[((i-1)*m+1):(i*m),1:m] <- kroneckerTtdTt[((i-1)*m+1):(i*m),1:m] + Tt

      prodTtP1<-model$P1
      dProdTtP1<-matrix(solve(kroneckerTtTt,kroneckerTtdTt%*%c(prodTtP1)),m,m)
      dV[1,i]<-dProdTtP1[1]
      for(t in 2:n){
        dProdTtP1<-Tt%*%dProdTtP1+dTt%*%prodTtP1
        dV[t,i]<-dProdTtP1[1]
        prodTtP1<-Tt%*%prodTtP1
      }
    }
  }
  if(q>0){
    Rt<-model$R
    dRtRt<-matrix(0,m,m)
    for(i in 1:q){
      dRtRt[]<-0
      dRtRt[i+1,]<-Rt
      dProdTtP1<-matrix(solve(kroneckerTtTt,c(dRtRt+t(dRtRt))),m,m)
      dV[1,p+i]<-dProdTtP1[1]
      for(t in 2:n){
        dProdTtP1<-Tt%*%dProdTtP1
        dV[t,p+i]<-dProdTtP1[1]
      }
    }
  }
  dV
}

