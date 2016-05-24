##############################################################
#' @name kobeTrks
#' @title Creates tracks for plotting on phase plot
#' @description
#' 
#' Creates tracks for plotting in the Kobe Phase Plot
#' 
#' @aliases kobeTrks-method kobeTrks,numeric,numeric-method kobeTrks,data.frame,missing-method
#' 
#'
#' @param   stock   an object of class \code{vector,data.frame,...}
#' @param   harvest an object of class \code{vector,data.frame,...}
#' @param   prob    a numeric vector with probabilities
#' @param   na.rm   logical
#' @export
#' @docType methods
#' @rdname  kobeTrks-method
#'
#' @examples
#' \dontrun{kobeTrks()}
#'
setMethod('kobeTrks', signature(stock="numeric",harvest="numeric"),
      function(stock,harvest,prob=c(0.25,0.5,0.75),na.rm=FALSE){
                                   
       res=kobeTrksFn(data.frame(stock=stock,harvest=harvest),prob=prob,na.rm=na.rm)
        
     return(res)})
 
setMethod('kobeTrks', signature(stock='data.frame',harvest="missing"),
     function(stock,harvest="missing",prob=c(0.25,0.5,0.75),na.rm=FALSE){
     res=kobeTrksFn(stock,prob=prob,na.rm=na.rm)
       
     return(res)})

kobeTrksFn=function(x,prob=c(0.25,0.5,0.75),na.rm=FALSE){
  res=apply(x,2,quantile,prob=prob,na.rm=na.rm)
  if (length(prob)>1)
     res=data.frame("Percentile"=dimnames(res)[[1]],res)
  else
     res=data.frame("Percentile"=paste(prob*100,"%",sep=""),t(res))
  
  res}
