#' @title Standardization of the qualitative variables
#' 
#' @param X.quali : a factor or a data frame with several factors
#' @param metric : the metric to be used, i.e. each category is weighted by the inverse of the square-root of its relative frequency
#'
#' @return Xdisj.sd : a standardized matrix with as many columns as categories associated with the qualitative variables.
#' 
#' @export
#' 
 
stand_quali =  function (X.quali, metric="chisq") 
{
  
  if (is.numeric(X.quali))  stop("All variables in X.quali must be factors")
  Xdisj<-tabdisj(X.quali)
  D<-apply(Xdisj,2,mean)
  Xdisj.sd<-Xdisj%*%diag(sqrt(1/D))
  colnames(Xdisj.sd)<-colnames(Xdisj)
  return(Xdisj.sd)
  
}
