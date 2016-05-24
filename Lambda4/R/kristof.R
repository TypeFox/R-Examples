#'  Compute Kristof Coefficient
#'
#' @description A reliability coefficient used for tests that are easily split into three parts.
#'
#' @return
#' \item{kristof}{The Kristof estimate of reliability.}
#' \item{Split}{The split used to obtain the reliability estimate.}
#'
#'
#' @param x Can be either a data matrix or a covariance matrix
#' @param split.method Specify method for splitting items?
#' @param missing How to handle missing values.
#' @param standardize When TRUE Results are standardized by using the correlation matrix instead of the covariance matrix for computation.
#'
#' @author Tyler Hunt \email{tyler@@psychoanalytix.com}
#' @references
#' Kristof, W. (1974). Estimation of reliability and true score variance from a split of a test into three arbitrary parts. Psychometrika, 39(4), 491-499.
#' @examples
#' kristof(Rosenberg, split.method="triplet")
#' @export
kristof<-function(x, split.method="triplet", missing="complete", standardize=FALSE){

 n <- dim(x)[1]
 p <- dim(x)[2]
	
 sigma <- impute.cov(x, missing)
 
 if(standardize==TRUE){
   sigma <- cov2cor(sigma)
 }
 
  if(split.method[1]=="triplet") {
  	rr<-rep(c(-1,0,1),ceiling(p/3))[1:p]
  	t1<-matrix((abs(rr)-1)*-1, ncol=1)
 	t2<-matrix((abs(rr/2)+rr/2), ncol=1)
 	t3<-matrix((abs(rr/2)-rr/2), ncol=1)
 }
 
 if(split.method[1]=="random") {
 	rr<-round(runif(p, -1.5, 1.5))
 	t1<-matrix((abs(rr)-1)*-1, ncol=1)
 	t2<-matrix((abs(rr/2)+rr/2), ncol=1)
 	t3<-matrix((abs(rr/2)-rr/2), ncol=1)
 }
 
 if(split.method[1]=="evenly.random") {
 	rr<-sample(rep(c(-1,0,1),ceiling(p/3))[1:p])
  	t1<-matrix((abs(rr)-1)*-1, ncol=1)
 	t2<-matrix((abs(rr/2)+rr/2), ncol=1)
 	t3<-matrix((abs(rr/2)-rr/2), ncol=1)
 }
 
  if(split.method[1]==-1 | split.method[1]==0 | split.method[1]==1) 
  	rr<-split.method 
  	t1<-matrix((abs(rr)-1)*-1, ncol=1)
 	t2<-matrix((abs(rr/2)+rr/2), ncol=1)
 	t3<-matrix((abs(rr/2)-rr/2), ncol=1)
  	
 
 t1t<-t(t1)
 t2t<-t(t2)
 t3t<-t(t3)
 
kristof<-(t1t%*%sigma%*%t2*t1t%*%sigma%*%t3+
t1t%*%sigma%*%t2*t2t%*%sigma%*%t3+
t1t%*%sigma%*%t3*t2t%*%sigma%*%t3)^2/(
t1t%*%sigma%*%t2*t1t%*%sigma%*%t3*t2t%*%sigma%*%t3*sum(sigma))

Split=cbind(t1,t2,t3)

result<-list(kristof=kristof, Split=Split)

class(result)="kristof"

return(result)

}
