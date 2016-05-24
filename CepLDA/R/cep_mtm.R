#' @name cep.mtm
#' @aliases cep.mtm
#' @title Multitaper Estimation of Cepstral Coefficients and the Log-Spectrum
#' 
#' @description
#' Returns multitaper estimated cepstra coefficients and log-spectrum for univariate time series.
#' 
#' @usage
#' cep.mtm(x,nw,k)
#' 
#' @param x Univariate time series of length N.
#' @param nw Width of tapers. Default is set to 4
#' @param k Number of tapers. Default is set to 7
#' @export
#' @importFrom stats fft
#' @import stats
#' @return a list with 4 elements
#' \item{quef}{Quefencies.}
#' \item{cep}{Raw cepstra coefficients from 0:(N-1)}
#' \item{freq}{Frequencies between 0 and 1}
#' \item{lspec}{Log-spectrum.}
#' @author Robert Krafty \email{<rkrafty@@pitt.edu>}
#' @examples
#' ## simulate a time series
#' N = 500 #length of each series
#' dat <- r.cond.ar2(N=N,nj=1,r.phi1=c(.01,.7),r.phi2=c(-.12,-.06),r.sig2=c(.3,3))$X
#' 
#' ## Fit multiaper
#' cep <- cep.mtm(dat)
#' 
#' ## Plot the cepstral coefficients
#' plot(cep$quef, cep$cep)
#' 
#' ## Plot the log spectrum
#' plot(cep$freq, cep$lspec, type="l")

cep.mtm <-function(x, nw=4, k=7){
        n = length(x)
        lpa = lperd.mtm(x, nw, k)
        lp1 = lpa$lspec
        cp = fft(c(lp1[n],lp1[-n]))*sqrt(2)/n
        cp[1] = cp[1]/sqrt(2)
        z=list(quef=0:(n-1), cep=Re(cp), freq=lpa$freq, lspec=lpa$lspec)
        z
}