#' Compute Raju Coefficient
#'
#' @param x Can be either a data matrix or a covariance matrix
#' @param split.method Specify method for splitting items.
#' @param missing How to handle missing values.
#' @param standardize When TRUE Results are standardized by using the correlation matrix instead of the covariance matrix for computation.
#'
#' @author Tyler Hunt \email{tyler@@psychoanalytix.com}
#' @examples
#' raju(Rosenberg, split.method="even.odd")
#' @export
raju<-function(x, split.method="even.odd", missing="complete", standardize=FALSE){

 n <- dim(x)[1]
 p <- dim(x)[2]
	
 sigma <- impute.cov(x, missing)
 
 if(standardize==TRUE){
   sigma <- cov2cor(sigma)
 }
 
 if(split.method[1]=="even.odd") 
 	t1t.split<-rep(c(1,0),ceiling(p/2))[1:p] 
 if(split.method[1]=="random") 
 	t1t.split<-round(runif(p))
 if(split.method[1]=="evenly.random") 	
 	t1t.split<-sample(rep(c(1,0),ceiling(p/2))[1:p])
 if(split.method[1]==1 | split.method[1]==0) 
 	t1t.split<-split.method 
 if(length(t1t.split)!=p)
 	warning("The length of split is not the same as the number of items")

 Split<-t1t.split
 t1<-matrix(Split, ncol=1)
 t1t<-t(t1)
 t2<-(t1-1)*-1
 t2t<-t(t2)
 
 lisq<-(sum(t1)/length(t1))^2+(sum(t2)/length(t2))^2
 
 raju.est<-(sum(sigma)-(t1t%*%sigma%*%t1+t2t%*%sigma%*%t2))/((1-lisq)*sum(sigma))

 result<-c(raju.est=raju.est)
 class(result)<-c("raju")
 return(result)

}

