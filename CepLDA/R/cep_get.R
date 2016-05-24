#' @name cep.get
#' @aliases cep.get
#' @title Obtain Data Frame to be Used in \code{Lopt.get} and \code{predict.ceplda}
#' 
#' @description
#' Returns a data frame containing raw cepstra coefficients and the group membership from multiple time seres.
#' 
#' @usage
#' cep.get(y,x,nw,k)
#' 
#' @param y n-vector indicating  group membership
#' @param x \emph{N} by \emph{n} matrix containing \emph{n} time series, each with length \emph{N}. 
#' @param nw Width of tapers used in multitaper spectral estimation. Default is set to 4
#' @param k Number of tapers used in multitaper spectral estimation. Default is set to 7
#' @export
#' @return 
#' \item{D.hat}{Data frame containing group information and raw cepstral coefficients.}
#' @references Krafty, RT(2016) Discriminant Analysis of Time Series in the Presence of Within-Group Spectral Variability.  \emph{Journal of Time series analysis} 
#' @author Zeda Li \email{<zeda.li@@temple.edu>}
#' @seealso 
#'  \code{\link{predict.ceplda}}, \code{\link{Lopt.get}}
#' @examples 
#' ## Simulate dataset
#' nj = 50  #number of series in training data
#' N = 500  #length of time series
#' data1 <- r.cond.ar2(N=N,nj=nj,r.phi1=c(.01,.7),r.phi2=c(-.12,-.06),r.sig2=c(.3,3))
#' data2 <- r.cond.ar2(N=N,nj=nj,r.phi1=c(.5,1.2),r.phi2=c(-.36,-.25),r.sig2=c(.3,3))
#' data3 <- r.cond.ar2(N=N,nj=nj,r.phi1=c(.9,1.5),r.phi2=c(-.56,-.75),r.sig2=c(.3,3))
#' dat <- cbind(data1$X,data2$X,data3$X)
#' y <- c(rep(1,nj),rep(2,nj),rep(3,nj))
#' data.cep <- cep.get(y,dat,4,7)
#' dim(data.cep)

cep.get <- function(
        y,              
        x,              
        nw = 4,        
        k = 7           
){
        if(dim(x)[2] != length(y))
                stop("\n Number of time series and group information (y) must be the same \n")
        #################################################
        ## get random spectra and cepstral coefficeints
        n <- dim(x)[2]
        N <- dim(x)[1]
        log.spec.hat <- matrix(0,nrow=n,ncol=N)
        cep.hat <- matrix(0,nrow=n,ncol=N)
        for(j in 1:n){
                temp = cep.mtm(x[,j],nw=nw,k=k)
                log.spec.hat[j,] = temp$lspec
                cep.hat[j,]= temp$cep
        }
        D.hat <- cbind(data.frame(cep.hat),y)
        colnames(D.hat)[1:N]= paste("C", 0:(N-1), sep="")
        return(D.hat)
}