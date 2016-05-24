#' @name predict.ceplda
#' @aliases predict.ceplda
#' @title Classify Multivariate Time Series
#' 
#' @description
#' Classify time series.  Run as part of \code{cep.lda}, and can be run seperatly after running \code{cep.lda}.
#' @param object Object of class "ceplda".
#' @param newdata Data frame of cases to be classified. Data frame can be obtained by \code{cep.get}.
#' @param ... argument based from or to other methods.
#' @details This function is a method for the generic function \code{predict()} for class "ceplda". 
#' It can be invoked by calling \code{predict(x)} for an object x of the appropriate class, 
#' or directly by calling \code{predict.lda(x)} regardless of the class of the object. Details usage of this function
#' is similar to \code{predict.lda(MASS)}.
#' @method predict ceplda
#' @export
#' @importFrom stats predict
#' @return 
#' List with components
#' \item{class}{Classification result (a factor).}
#' \item{posterior}{Posterior class probabilities.}
#' \item{x}{Scores from test data.}
#' @seealso
#'  \code{\link{cep.lda}}, \code{\link{predict.lda}}
#' @examples 
#' ## See cep.lda for predicting new data simultaneously while fitting a model to training data.
#' ## Below is predicting new data after fitting a model to the training data.  
#' 
#' ## Simulate training data
#' nj = 50  #number of series in training data
#' N = 500  #length of time series
#' traindata1 <- r.cond.ar2(N=N,nj=nj,r.phi1=c(.01,.7),r.phi2=c(-.12,-.06),r.sig2=c(.3,3))
#' traindata2 <- r.cond.ar2(N=N,nj=nj,r.phi1=c(.5,1.2),r.phi2=c(-.36,-.25),r.sig2=c(.3,3))
#' traindata3 <- r.cond.ar2(N=N,nj=nj,r.phi1=c(.9,1.5),r.phi2=c(-.56,-.75),r.sig2=c(.3,3))
#' train <- cbind(traindata1$X,traindata2$X,traindata3$X)
#' y <- c(rep(1,nj),rep(2,nj),rep(3,nj))
#' 
#' ## Fit the discriminant analysis
#' fit <- cep.lda(y,train)
#' 
#' ## Simulate test data
#' testdata1 <- r.cond.ar2(N=N,nj=nj,r.phi1=c(.01,.7),r.phi2=c(-.12,-.06),r.sig2=c(.3,3))
#' testdata2 <- r.cond.ar2(N=N,nj=nj,r.phi1=c(.5,1.2),r.phi2=c(-.36,-.25),r.sig2=c(.3,3))
#' testdata3 <- r.cond.ar2(N=N,nj=nj,r.phi1=c(.9,1.5),r.phi2=c(-.56,-.75),r.sig2=c(.3,3))
#' test <- cbind(testdata1$X,testdata2$X,testdata3$X)
#' yTest <- c(rep(1,nj),rep(2,nj),rep(3,nj))
#' 
#' 
#' ## Classifty new data
#' pre <- predict(fit,cep.get(yTest,test))
#' mean(pre$class == yTest)
#' table(yTest,pre$class)



predict.ceplda <- function(object, newdata,...){
  if(!inherits(object,"ceplda"))
    stop("Object must be of class 'ceplda'")
  if(length(object$C.lda)!=10)
    stop("Object muse be from lda without cross-validation")
  return(predict(object$C.lda, newdata,...))
}