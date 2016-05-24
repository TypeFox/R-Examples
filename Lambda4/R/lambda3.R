#' Compute Guttman's Lambda 3 Coefficient (Coefficent Alpha)
#'
#' @description Often recognized as Cronbach's alpha, Guttman's Lambda 3 can be used to estimate reliability when the data can be split in parallel forms.
#' 
#' @return
#' \item{lambda3}{The unstandardized and standardized lambda3 estimate.}
#' \item{item.stats}{If the input data was a covariance matrix then this is a table of reliability estimates if an item was dropped.  If the input data is a data frame then the mean, standard deviation, and number of observations are also included.}
#' \item{items}{The number of items.}
#' \item{item.stats.max}{The maximum number of item to display the item.stats table (user specified).}
#' 
#' @param x Can be either a data matrix or a covariance matrix
#' @param item.stats.max items statistics shown if the number of items are less than this value.
#' @param missing how to handle missing values.
#' 
#' @references
#' Cronbach L (1951). "Coefficient Alpha and the Internal Structure of Tests." Psychometrika, 16, 297-334.
#' Guttman L (1945). "A Basis for Analyzing Test-Retest Reliability." Psychometrika, 10, 255-282.
#'
#' @author Tyler Hunt \email{tyler@@psychoanalytix.com}
#' @examples 
#' lambda3(Rosenberg)
#' @export


lambda3<-function(x, item.stats.max=12, missing="complete"){
	
  n <- dim(x)[1]
  p <- dim(x)[2]
    
  sigma <- impute.cov(x, missing)
 
  sigma.cor<-cov2cor(sigma)
 
  Obs<-colSums(!is.na(x))
  Mean<-colMeans(x, na.rm=TRUE)
  SD<-sapply(x, sd, na.rm=TRUE)
 
  onerow<-rep(1,p)
  onerow<-t(onerow)
  onevector<-t(onerow)
 
  Unstandardized<-(p/(p-1))*(1-(onerow%*%diag(sigma)/(onerow%*%sigma%*%onevector)))
  Standardized<-(p/(p-1))*(1-(onerow%*%diag(sigma.cor)/(onerow%*%sigma.cor%*%onevector)))
  items<-p
  lambda3<-data.frame(Unstandardized, Standardized)
 
  If.Dropped<-rep(NA,p)
  for(i in 1:p){
    onerow.d<-rep(1,(p-1))
    onerow.d<-t(onerow.d)
    onevector.d<-t(onerow.d)
 	  sigma.d<-sigma[-i,-i]
 	  If.Dropped[i]<-(p/(p-1))*(1-(onerow.d%*%diag(sigma.d)/(onerow.d%*%sigma.d%*%onevector.d)))
  }
  
  if(n != p) {
    item.stats<-data.frame(Mean,SD,Obs,If.Dropped, row.names=(colnames(x)))
  }
  else{
    item.stats=data.frame(If.Dropped, row.names=(colnames(x)))
  }
  rownames(item.stats)=colnames(x)
    
  result<-list(lambda3=lambda3, item.stats=item.stats, items=items, item.stats.max=item.stats.max) 

 
  class(result) <- c("lambda3")
  return(result)
}
