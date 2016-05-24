#' Compute Quantile Lambda 4
#'
#'
#' @description Quantile maximize lambda4 is a statistic that can be used in most measurement situations.  In particular this function generates a vector t of length equal to the number of items.  Each value in the vector consists of either a +1 or -1 (randomly generated).  Next, in a random order each value in the t-vector is switched.  The value kept (+1 or -1) is the value that resulted in the highest reliability estimate.  This procedure is repeated by default 1000 times but can also be user specified.  The user can then specify the quantile of this vector but it defaults to .5.
#'
#'
#' @param x Can be either a data matrix or a covariance matrix
#' @param starts How many split-half reliability estimates used
#' @param quantiles The quantiles of the generated splits.  It defaults to .5 because it makes the most sense at this time.  (The simulation manuscript is under review).
#' @param missing How to handle missing values.
#' @param show.lambda4s If TRUE then Shows the vector of lambda4s if FALSE then the vector is hidden
#' @param standardize Results are standardized by using the correlation matrix instead of the covariance matrix for computation.
#'
#'
#' @references
#' Guttman L (1945). "A Basis for Analyzing Test-Retest Reliability." Psychometrika, 10, 255-282. 
#' 
#' Callender J, Osburn H (1977). "A Method for Maximizing and Cross-Validating Split-Half Reliability Coefficients." Educational and Psychological Measurement, 37, 819-826.
#' 
#' Callender J, Osburn H (1979). "An Empirical Comparison of Coefficient Alpha, Guttman's Lambda2 and Msplit Maximized Split-Half Reliability Estimates." Journal of Educational Measurement, 16, 89-99. 
#' Guttman L (1945). "A Basis for Analyzing Test-Retest Reliability." Psychometrika, 10, 255-282.
#' 
#' Callender J, Osburn H (1977). "A Method for Maximizing and Cross-Validating Split-Half Reliability Coefficients." Educational and Psychological Measurement, 37, 819-826. 
#' 
#' Callender J, Osburn H (1979). "An Empirical Comparison of Coefficient Alpha, Guttman's Lambda2 and Msplit Maximized Split-Half Reliability Estimates." Journal of Educational Measurement, 16, 89-99.
#' 
#' Sijtsma K (2009). "On the Use, Misuse, and Very Limited Usefulness of Cronbach's Alpha." Psychometrika, 74(1), 107-120.
#'
#' @return
#' \item{lambda4.quantile}{The user specified quantile value of the vector of maximized split-reliability}
#' \item{lambda4.optimal}{Maximum split-half reliability (Maximized Lambda4}
#' \item{l4.vect}{A vector of lambda4 (split-half reliability) calculations}
#'
#' @examples
#' quant.lambda4(Rosenberg, starts=1000, quantile=c(.05,.5,.95))
#'
#' @author Tyler Hunt \email{tyler@@psychoanalytix.com}
#' 
#' @export


quant.lambda4<-function(x, starts=1000, quantiles=.5, missing="complete", show.lambda4s=FALSE, standardize=FALSE){

  l4.vect<-rep(NA, starts)

#Determines if x is a covariance or data matrix and establishes a covariance matrix for estimation.

  sigma <- impute.cov(x, missing)

  if(standardize==TRUE){
    sigma <- cov2cor(sigma)
  }

  items<-ncol(sigma)

#Creates an empty matrix for the minimized tvectors
  splitmtrx<-matrix(NA, nrow=items, ncol=starts)

# creates the row and column vectors of 1s for the lambda4 equation.
  onerow<-rep(1,items)
  onerow<-t(onerow)
  onevector<-t(onerow)
  f<-rep(NA,starts)

  sigma0<-sigma
  diag(sigma0)<-0

  for(y in 1:starts){

#Random number generator for the t-vectors
    trow<-(round(runif(items, min=0, max=1))-.5)*2
    trow<-t(trow)
    tvector<-t(trow)

#Creating t vector and row
    tk1<-(tvector)
    tk1t<-t(tk1)
    tk2<-(trow)
    tk2t<-t(tk2)

#Decision rule that determines which split each item should be on.  Thus minimizing the numerator.

    random.order<-sample(1:items)
    for (o in 1:items){
      oi<-sigma0[,random.order[o]]
      fi<-oi%*%tk1
      if (fi <  0) {tk1[random.order[o],1]<-  1}
      if (fi >= 0) {tk1[random.order[o],1]<- -1}
    }

    t1<-(1/2)*(tk1+1)
    fk1<-tk1t%*%sigma0%*%tk1
    t1t<-t(t1)
    t2<-(1-t1)
    t2t<-t(t2)

    f[y]=fk1
    splitmtrx[,y]<-t1

    l4.vect[y]<-(4*(t1t%*%sigma%*%t2))/(onerow%*%sigma%*%onevector)
  }
  quants<-quantile(l4.vect, quantiles)

  lambda4.quantile=quants

  result<-list(lambda4.quantile=lambda4.quantile,
              lambda4s=l4.vect, 
              show.lambda4s=show.lambda4s,
              quantiles=quantiles)
  
  class(result)<-c("quant.lambda4")

  return(result)
}
