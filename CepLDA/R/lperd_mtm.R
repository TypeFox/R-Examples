#' @name lperd.mtm
#' @aliases lperd.mtm
#' @title Multitaper Log-Periodogram
#' @description
#' Obtain multivariate log-periodogram.  Not often used by itself, but is used as part of \code{cep.mtm}.
#' 
#' @usage
#' lperd.mtm(x,nw,k)
#' 
#' @param x \emph{N} by \emph{n} matrix containing \emph{n} training time series each with length \emph{N}.
#' @param nw Width of tapers used in multitaper spectral estimation. 
#' @param k Number of tapers used in multitaper spectral estimation. 
#' @export
#' @importFrom stats ts
#' @return 
#' \item{freq}{Fourier frequencies from 0 to 1.}
#' \item{lspec}{Log-periodogram.}
#' \item{spec}{Periodigram on the natural scale.}
#' @references Krafty, RT (2016) Discriminant Analysis of Time Series in the Presence of Within-Group Spectral Variability.  \emph{Journal of Time series analysis} 
#' @author Robert Krafty \email{rkrafty@pitt.edu}
#' 
#' @seealso
#'  \code{\link{cep.mtm}}, \code{\link{cep.lda}} 


lperd.mtm <-function(x, nw, k){
        ###### Univarite log-periodograms
        ## Get (1) spec, (2) log-spec, 
        ##    
        mtm = spec.mtm(ts(x), nw=nw, k=k, nFFT=length(x), returnZeroFreq=FALSE, plot=FALSE)
        #usehalf = (1:length(mtm$freq)) %% 2 == 0
        perdhalf = mtm$spec #[usehalf]
        lperdhalf = log(perdhalf)
        fr = mtm$freq #[usehalf]
        frq = (1:length(x))/length(x)
        ####### get the bias adjusted log periodograms
        if (fr[length(fr)] == .5) {
                lspec = c(lperdhalf[-length(perdhalf)],
                          lperdhalf[length(perdhalf)],
                          rev(lperdhalf[-length(perdhalf)]),
                          log(mean(perdhalf[1:5])))
                spec = c(perdhalf[-length(perdhalf)],
                         perdhalf[length(perdhalf)],
                         rev(perdhalf[-length(perdhalf)]),
                         mean(perdhalf[1:5]))
        }
        else {
                lspec = c(lperdhalf,  rev(lperdhalf), 
                          log(mean(perdhalf[1:5])))
                spec = c(perdhalf,  rev(perdhalf), mean(perdhalf[1:5]) )
        }
        ####### finish
        z=list(freq=frq, lspec=lspec, spec=spec)
}
