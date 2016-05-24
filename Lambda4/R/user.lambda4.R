#' Compute User Specified Lambda 4 (Split-Half)
#' 
#' @param x Can be either a data frame or a covariance matrix.
#' @param split.method Specify method for splitting items.
#' @param item.stats If TRUE then item statistics are provided in the output.
#' @param missing How to handle missing values.
#' 
#' @references
#' Guttman L (1945). "A Basis for Analyzing Test-Retest Reliability." Psychometrika, 10, 255-282.
#' @author Tyler Hunt \email{tyler@@psychoanalytix.com}
#' 
#' @examples
#' user.lambda4(Rosenberg)
#' user.lambda4(Rosenberg, c(0, 1, 1, 0, 1, 1, 0, 1, 0, 0))
#' 
#' @export


user.lambda4<-function(x, split.method="even.odd", item.stats=FALSE, missing="complete"){

 #number of variables
 nvar<-dim(x)[2] 
 #number of participants
 n<-dim(x)[1]
 
 #Determines if x is a covariance or a data matrix and establishes a covariance matrix for estimation.
 p <- dim(x)[2]
 sigma <- impute.cov(x, missing)
 
 if(split.method[1]=="even.odd") t1t.split<-rep(c(1,0),ceiling(nvar/2))[1:nvar] 
 if(split.method[1]=="random") t1t.split<-round(runif(nvar))
 if(split.method[1]=="evenly.random") t1t.split<-sample(rep(c(1,0),ceiling(nvar/2))[1:nvar])
 if(split.method[1]==1 | split.method[1]==0) t1t.split<-split.method 
 if(length(t1t.split)!=nvar)
 	warning("The length of split is not the same as the number of items")
 Split<-t1t.split
 Obs<-colSums(!is.na(x))
 

  if(n != p & item.stats == TRUE){
 	  Mean<-round(colMeans(x, na.rm=TRUE),digits=2)
 	  SD<-round(sapply(x,sd, na.rm=TRUE), digits=2)
 	  Item.Statistics<-data.frame(Mean,SD,Obs, row.names=(colnames(x)))
	  }
  else{
   Item.Statistics<-NULL
    }
  if(n == p & item.stats == TRUE){
    warning("Item statistics cannot be provided from a covariance matrix", call.=FALSE)
    Item.Statistics=NULL
    }


 t1t.split<-t(t1t.split)
 t2.split<-(t(t1t.split)-1)*-1
 
 onerow<-rep(1, nvar)
 onerow<-t(onerow)
 onevector<-t(onerow)

 lambda4<-(4*(t1t.split%*%sigma%*%t2.split))/sum(sigma)
 
 
 result<-list(lambda4=lambda4, Item.Statistics=Item.Statistics, Split=Split)
 
 class(result)=c("user.lambda4")
 return(result)
 }

