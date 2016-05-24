##############################################################
#' @name kobeSmry
#' @title Creates Summaries i.e. probability of being in the green quadrant
#' @description
#' 
#' Creates tracks for plotting in the Kobe Phase Plot
#' 
#' @aliases kobeSmry-method kobeSmry,numeric,numeric-method kobeSmry,data.frame,missing-method
#' 
#'
#' @param   stock   an object of class \code{vector,data.frame,...}
#' @param   harvest an object of class \code{vector,data.frame,...}
#' @param   onlyColours \code{logical} to report only quadrant colours
#' @export
#' @docType methods
#' @rdname  kobeSmry-method
#'
#' @examples
#' \dontrun{kobeSmry()}
#'
setMethod('kobeSmry', signature(stock="numeric",harvest="numeric"),
    function(stock,harvest,onlyColours=FALSE){
                                  
    res=kobeSmryFn(data.frame(stock=stock,harvest=harvest),onlyColours=onlyColours)
      
    return(res)})

setMethod('kobeSmry', signature(stock='data.frame',harvest="missing"),
    function(stock,onlyColours=FALSE){
     
    res=kobeSmryFn(stock,onlyColours=onlyColours)
      
    return(res)})

kobeSmryFn=function(x,onlyColours=onlyColours){
  
  res =cbind(x,   kobeP(x$stock,x$harvest))
  res = with(res, data.frame(stock       =median(stock,       na.rm=TRUE),
                             harvest     =median(harvest,     na.rm=TRUE),
                             red         =mean(  red,         na.rm=TRUE),
                             yellow      =mean(  yellow,      na.rm=TRUE),
                             green       =mean(  green,       na.rm=TRUE),
                             overFished  =mean(  overFished,  na.rm=TRUE),
                             overFishing =mean(  overFishing, na.rm=TRUE),
                             underFished =1-mean(overFished,  na.rm=TRUE),
                             underFishing=1-mean(overFishing, na.rm=TRUE)))
  
  if (onlyColours) res=res[,3:5]
  return(res)}


