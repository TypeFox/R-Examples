#' @name r.cond.ar2
#' @aliases r.cond.ar2
#' @title Generate Random AR(2) Time Series
#' 
#' @description
#' Simulates multiple AR(2) time series.
#' 
#' @usage
#' r.cond.ar2(N,nj,r.phi1,r.phi2,r.sig2)
#' 
#' @param N  Length of the series
#' @param nj Number of series generated.
#' @param r.phi1 Range of first AR order coefficient. It is a vector contains minimum and maximum possible coefficients.
#' @param r.phi2 Range of second AR order coefficient. It is a vector contains minimum and maximum possible coefficients.
#' @param r.sig2 Range of conditional innovation variances. It is a vector contains minimum and maximum possible variances.
#' @export
#' @importFrom stats ts
#' @importFrom stats arima.sim
#' @importFrom stats runif
#' @return a list with 2 elements
#' \item{X}{\emph{N} by \emph{nj} matrix of time series}
#' \item{cep}{3 by \emph{nj} matrix of parameters (phi1, phi2, sig2)}
#' @references Krafty, RT (2016) Discriminant Analysis of Time Series in the Presence of Within-Group Spectral Variability.  \emph{Journal of Time series analysis} 
#' @author Robert Krafty \email{<rkrafty@@pitt.edu>}
#' @seealso
#'  \code{\link{cep.lda}} 
#' @examples
#' ## Simulate data
#' nj = 50  #number of series in training data
#' N = 500  #length of time series
#' data1 <- r.cond.ar2(N=N,nj=nj,r.phi1=c(.01,.7),r.phi2=c(-.12,-.06),r.sig2=c(.3,3))
#' data2 <- r.cond.ar2(N=N,nj=nj,r.phi1=c(.5,1.2),r.phi2=c(-.36,-.25),r.sig2=c(.3,3))
#' data3 <- r.cond.ar2(N=N,nj=nj,r.phi1=c(.9,1.5),r.phi2=c(-.56,-.75),r.sig2=c(.3,3))
#' data <- cbind(data1$X,data2$X,data3$X)
#' y <- c(rep(1,nj),rep(2,nj),rep(3,nj))  

r.cond.ar2 <- function(N, nj=1, r.phi1, r.phi2, r.sig2){
        # the notation here is consistant with the final paper
        # N is the length of a series
        # nj is the number of series
        # X is N x nj matrix of time series
        # parms is 3 x nj matrix of paramters (phi1, phi2, sig2)
        X = ts(matrix(0,ncol=nj,nrow=N))
        parms = matrix(0,ncol=nj,nrow=3)
        for(j in 1:nj){
                phi1 = runif(1,min=min(r.phi1), max=max(r.phi1) )
                phi2 = runif(1,min=min(r.phi2), max=max(r.phi2[2]) )
                sigma2 = runif(1,min=min(r.sig2), max=max(r.sig2[2]) )
                Xj = arima.sim(list(order=c(2,0,0),ar=c(phi1,phi2)),sd=sqrt(sigma2),n=N)
                X[,j] = Xj
                parms[,j] = c(phi1,phi2,sigma2)
        }
        z=list(X=X, parms=parms)
}



