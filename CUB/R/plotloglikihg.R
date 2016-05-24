#' @title Plot of the log-likelihood function of the IHG distribution
#' @aliases plotloglikihg
#' @description Plot the log-likelihood function of an IHG model fitted to a given absolute frequency distribution,
#'  over the whole support of the preference parameter. It returns also its argmax.
#' @usage plotloglikihg(m, freq)
#' @export plotloglikihg
#' @param m Number of ordinal categories
#' @param freq Vector of the absolute frequency distribution
#' @seealso \code{\link{loglikIHG}}
#' @examples
#' m<-7
#' freq<-c(828, 275, 202, 178, 143, 110, 101)
#' max<-plotloglikihg(m, freq)

plotloglikihg <-
function(m,freq){
  np<-1000
  ordinate<-rep(NA,np)
  ini<-1; fin<-np;
  thetavec<-(1:np)/np
  for(j in 1:np){
    ordinate[j]<-loglikihg(m,freq,thetavec[j])
  }
  plot(thetavec[ini:fin],ordinate[ini:fin],type="l",lwd=3,xlab=expression(theta),ylab="Log-likelihood function", 
       cex.main=0.9,main="Log-likelihood function for the Inverse HyperGeometric distribution")
  which.max(ordinate)/np ### MLE of theta
}
