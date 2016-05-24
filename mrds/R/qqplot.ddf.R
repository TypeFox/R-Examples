#' Q-Q plot, KS and CVM goodness of fit tests for distance detection functions
#'
#' Constructs a Q-Q plot for fitted model as a graphical picture of goodness of
#' fit and computes K-S and Cramer-VonMises goodness of fit tests for distance
#' sampling models based on single observer survey and double observer survey
#' with independent observer (io) and trial configurations.
#'
#' \code{pks} computes the p-value for the K-S test.  The function
#' \code{pcramer} was taken from the coda package.  It computes the p-value for
#' the CvM test.  Both \code{pks} and \code{pcramer} are used in
#' \code{qqplot.ddf} and need not be called by user. \code{qqplot.ddf} is
#' called from \code{ddf.gof} to evaluate model goodness of fit.
#'
#' @aliases qqplot.ddf pks pcramer
#' @usage qqplot.ddf(model, plot=TRUE, ...)
#'
#'        pcramer(q, eps = 1e-05)
#'
#'        pks(Dn,n)
#'
#' @param model fitted distance detection function model object
#' @param plot the Q-Q plot be plotted or just report statistics?
#' @param n sample size
#' @param Dn K-S statistic
#' @param q CvM statistic
#' @param eps small value that controls accuracy of p-value computation
#' @param \dots unspecified arguments passed to plot
#' @export
#' @return A list of goodness of fit related values: \item{edf}{matrix of lower
#'   and upper empirical distribution function values} \item{cdf}{fitted
#'   cumulative distribution function values} \item{ks}{list with K-S statistic
#'   (\code{Dn}) and p-value (\code{p})} \item{CvM}{list with CvM statistic
#'   (\code{W}) and p-value (\code{p})}
#' @author Jeff Laake
#' @seealso \code{\link{ddf.gof}}, \code{\link{cdf.ds}}
#' @references Burnham, K.P., S.T. Buckland, J.L. Laake, D.L. Borchers, T.A.
#'   Marques, J.R.B. Bishop, and L. Thomas. 2004.  Further topics in distance
#'   sampling. pp: 385-389. In: Advanced Distance Sampling, eds. S.T. Buckland,
#'   D.R.Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers, and L. Thomas.
#'   Oxford University Press.
#' @keywords utility
#' @importFrom graphics abline
qqplot.ddf <- function(model,plot=TRUE,...){

  fun <- function(x,z,lt){
    if(lt){
      length(z[z<x])
    }else{
      length(z[z<=x])
    }
  }

  if("ds" %in% class(model)){
    cdfvalues <- sort(cdf.ds(model)$fitted)
  }else if("io" %in% class(model) |
           "trial" %in% class(model) |
           "rem" %in% class(model) ){
    cdfvalues <- sort(cdf.ds(model$ds)$fitted)
  }else if("io.fi" %in% class(model)){
    data <- model$data
    data <- data[data$object %in% as.numeric(names(model$fitted)),]
    n <- length(data$distance)/2
    cdfvalues <- rep(0,n)
    for(i in 1:n){
      newdata <- data[data$object %in% as.numeric(names(model$fitted))[i],]
      cdfvalues[i] <- predict.io.fi(model,newdata=newdata,integrate=TRUE,
                                    int.range=newdata$distance[1])$fitted
      if(model$meta.data$left!=0){
        widthvalue <- predict.io.fi(model,newdata=newdata,integrate=TRUE,
                                    int.range=model$meta.data$width)$fitted
        cdfvalues[i] <- cdfvalues[i]/widthvalue
      }
    }
    if(model$meta.data$left==0){
      cdfvalues <- cdfvalues/predict.io.fi(model,integrate=TRUE)$fitted
    }
    cdfvalues <- sort(cdfvalues)
  }else if("trial.fi" %in% class(model)){
    data <- model$data
    data <- data[data$observer==1&data$object %in% 
                   as.numeric(names(model$fitted)),]
    n <- length(data$distance)
    cdfvalues <- rep(0,n)
    for(i in 1:n){
      newdata <- data[data$object %in% as.numeric(names(model$fitted))[i],]
      cdfvalues[i] <- predict.trial.fi(model,newdata=newdata,integrate=TRUE,
                                       int.range=newdata$distance[1])$fitted

      if(model$meta.data$left != 0){
        widthvalue <- predict.trial.fi(model,newdata=newdata,integrate=TRUE,
                                       int.range=model$meta.data$width)$fitted
        cdfvalues[i] <- cdfvalues[i]/widthvalue
      }
    }
    if(model$meta.data$left==0){
      cdfvalues <- cdfvalues/predict.trial.fi(model,newdata=data,
                                              integrate=TRUE)$fitted
    }
    cdfvalues <- sort(cdfvalues)
  }else if("rem.fi" %in% class(model)){
    data <- model$data
    data <- data[data$object %in% as.numeric(names(model$fitted)),]
    n <- length(data$distance)/2
    cdfvalues <- rep(0,n)
    for(i in 1:n){
      newdata <- data[data$object %in% as.numeric(names(model$fitted))[i],]
      cdfvalues[i] <- predict.rem.fi(model,newdata=newdata,integrate=TRUE,
                                     int.range=newdata$distance[1])$fitted

      if(model$meta.data$left!=0){
        widthvalue <- predict.rem.fi(model,newdata=newdata,integrate=TRUE,
                                  int.range=model$meta.data$width)$fitted
        cdfvalues[i] <- cdfvalues[i]/widthvalue
      }
    }
    if(model$meta.data$left==0){
      cdfvalues <- cdfvalues/predict.rem.fi(model,newdata=data,
                                            integrate=TRUE)$fitted
    }
    cdfvalues <- sort(cdfvalues)
  }

  n <- length(cdfvalues)
  lower.edf <- (unlist(sapply(cdfvalues,fun,z=cdfvalues,lt=TRUE)))/n
  upper.edf <- (unlist(sapply(cdfvalues,fun,z=cdfvalues,lt=FALSE)))/n

  if(plot){
    plot(upper.edf,cdfvalues,xlab="Empirical cdf",ylab="Fitted cdf",
         xlim=c(0,1),ylim=c(0,1),...)
    abline(0,1,...)
  }

  Dn <- max(max(abs(lower.edf-cdfvalues)),max(abs(upper.edf-cdfvalues)))
  W <- 1/(12*n) + sum((cdfvalues - ((1:n)-.5)/n)^2)
  return(list(edf=cbind(lower.edf,upper.edf),
              cdf=cdfvalues,
              ks=list(Dn=Dn,p=pks(Dn,n)),
              CvM=list(W=W,p=1-pcramer(W))
             )
        )
}
