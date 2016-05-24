#' Compute Angoff Coefficient
#'
#' @description Angoff's coefficient is most appropriately used for estimating reliability in tests that can be split into two parts with unequal lengths.  The calculation corrects for the inequality of length in the splits.  Angoff's coefficient is also believed to handle congeneric test structures relatively well.
#'
#'
#' @param x Can be either a data matrix or a covariance matrix
#' @param split.method Specify method for splitting items?
#' @param missing How to handle missing values.
#' @param standardize When TRUE Results are standardized by using the correlation matrix instead of the covariance matrix for computation.
#'
#' @return 
#' \item{angoff}{The estimate of reliability.}
#' \item{Split}{The split half key used to calculate angoff's coefficient.}
#' 
#' @author Tyler Hunt \email{tyler@@psychoanalytix.com}
#' 
#' @references
#' Feldt, L. S., & Charter, R. A. (2003). Estimating the reliability of a test split into two parts of equal or unequal length. Psychological Methods, 8(1), 102-109.
#' 
#' Sedere, M. U. And Feldt, L. S. (1977), The Sampling Distributions Of The Kristof Reliability Coefficient, The Feldt Coefficient, And Guttman's Lambda-2. Journal Of Educational Measurement, 14: 53-62. 
#' 
#' Feldt, L. S. (1975). Estimation of the reliability of a test divided into two parts of unequal length. Psychometrika, 40, 557-561.
#' 
#' Angoff, W. H. (1953). Test reliability and effective test length. Psychometrika, 18, 1-14.
#' 
#' @examples 
#' angoff(Rosenberg, split.method="even.odd", missing="complete", standardize=FALSE)
#' @export
angoff<-function(x, split.method="even.odd", missing="complete", standardize=FALSE){

 n <- dim(x)[1]
 p <- dim(x)[2]
	
 sigma <- impute.cov(x, missing)
 
 if(standardize==TRUE){
   sigma <- cov2cor(sigma)
 }
 
 if(split.method[1]=="even.odd") t1t.split<-rep(c(1,0),ceiling(p/2))[1:p] 
 if(split.method[1]=="random") t1t.split<-round(runif(p))
 if(split.method[1]=="evenly.random") t1t.split<-sample(rep(c(1,0),ceiling(p/2))[1:p])
 if(split.method[1]==1 | split.method[1]==0) t1t.split<-split.method 
 if(length(t1t.split)!=p)
 	warning("The length of split is not the same as the number of items")
 Split<-t1t.split
 
 t1<-matrix(Split, ncol=1)
 t1t<-t(t1)
 t2<-(t1-1)*-1
 t2t<-t(t2)
 
angoff<-4*(t1t%*%sigma%*%t2)/(sum(sigma)-((t1t%*%sigma%*%t1-t2t%*%sigma%*%t2)/sqrt(sum(sigma)))^2)

result<-list(angoff=angoff, Split=Split)
class(result)<-c("angoff")
return(result)

}
