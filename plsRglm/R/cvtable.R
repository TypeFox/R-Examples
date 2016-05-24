cvtable <- function(x,verbose=TRUE,...)
{
  if(class(x) %in% c("summary.cv.plsRmodel","summary.cv.plsRglmmodel")){
  if(class(x)=="summary.cv.plsRmodel"){return(cvtable.plsR(x,verbose=verbose))} 
  if(class(x)=="summary.cv.plsRglmmodel"){return(cvtable.plsRglm(x,verbose=verbose))}} else {
  stop("cvtable must be applied to a summary.cv.plsRmodel or summary.cv.plsRglmmodel object")
  }
}
