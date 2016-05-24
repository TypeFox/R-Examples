#' @name Lopt.get
#' @aliases Lopt.get
#' @title Optimal Choice of \eqn{L} 
#' @description
#' Data driven selection of the number of cepstral coefficients (\eqn{L}) via leave-one-out cross-validation.
#' 
#' @usage
#' Lopt.get(data,mcep)
#' 
#' @param data Data frame containing cepstral coefficients and group information. Obtained from \code{cep.get} or \code{cep.lda}.
#' @param mcep Maximum number of cepstral coefficient considerd. Default is set to 10
#' @export
#' @importFrom stats formula
#' @return 
#' \item{Lopt}{Optimal number of cepstral coefficients}
#' @references Krafty, RT (2016) Discriminant Analysis of Time Series in the Presence of Within-Group Spectral Variability.  \emph{Journal of Time series analysis} 
#' @author Robert Krafty \email{rkrafty@pitt.edu}
#' 
#' @seealso
#'  \code{\link{cep.lda}}, \code{\link{cep.get}} 
#' @examples
#' ## Simulate data
#' ntrain = 50  #number of series in training data
#' Nlength = 500 #length of each series
#' set.seed(2016)
#' traindata1 <- r.cond.ar2(N=Nlength,nj=ntrain,r.phi1=c(.01,.7),r.phi2=c(-.12,-.06),r.sig2=c(.3,3))
#' traindata2 <- r.cond.ar2(N=Nlength,nj=ntrain,r.phi1=c(.5,1.2),r.phi2=c(-.36,-.25),r.sig2=c(.3,3))
#' traindata3 <- r.cond.ar2(N=Nlength,nj=ntrain,r.phi1=c(.9,1.5),r.phi2=c(-.56,-.75),r.sig2=c(.3,3))
#' train <- cbind(traindata1$X,traindata2$X,traindata3$X)
#' ## group information
#' y <- c(rep(1,ntrain),rep(2,ntrain),rep(3,ntrain))
#' dat <- cep.get(y,train)
#' Lopt.get(dat,10)

Lopt.get <- function(data, mcep=10){
        cvK <- array(0,dim=mcep)
        for(k in 1:mcep){
                b <- as.formula(paste("y ~ ",paste(colnames(data[,1:(k+1)]), collapse="+"),sep = ""))
                C.lda.pred <- lda(b , data=data, CV=TRUE)
                cvK[k] <- mean(C.lda.pred$class==data$y)
        }
        Lopt <- min(which(cvK == max(cvK)))
        Lopt
}
