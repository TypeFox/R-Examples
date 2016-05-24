#' @import ggplot2 

utils::globalVariables(c("x","y","fill"))

setGeneric('kobePhase',  function(object,...)         standardGeneric('kobePhase'))

### provide a back drop on which to overlay data
kobePhaseFn=function(object,xlim,ylim){    
  quads<- rbind(data.frame(x=c(-Inf,-Inf,Inf,Inf), y=c(-Inf,Inf,Inf,-Inf), fill=as.factor("yellow")),
                data.frame(x=c(   1,   1,Inf,Inf), y=c(-Inf,  1,  1,-Inf), fill=as.factor("green")),
                data.frame(x=c(-Inf,-Inf,  1,  1), y=c(   1,Inf,Inf,   1), fill=as.factor("red")))
  
  p=ggplot(object)+geom_polygon(data=quads,aes(x,y,fill=fill)) +
    scale_fill_manual(values = c("yellow","green","red"), guide="none") +
    ylab(expression(F/F[MSY]))        +
    xlab(expression(SSB/B[MSY]))      +
    scale_y_continuous(limits=ylim)   +
    scale_x_continuous(limits=xlim)
  
  invisible(p)}

##############################################################
#' kobePhase
#'
#' @name kobePhase
#' 
#' @description
#' Creates The Kobe Phase Plot
#' @aliases kobePhase-method kobePhase,missing-method kobePhase,data.frame-method
#'
#' @param   object an object of class \code{missing,data.frame,...}
#' @param   xlim a numeric vector with x-axis limits, by default is c(0.2) 
#' @param   ylim a numeric vector with y-axis limits, by default is the same as xlim 
#' @export
#' @docType methods
#' @rdname  kobePhase-method
#'
#' @examples
#' \dontrun{kobePhase()}
setMethod('kobePhase', signature(object='missing'),
  function(object,xlim=c(0,2),ylim=xlim){
    
       invisible(kobePhaseFn(NULL,xlim,ylim))})

setMethod('kobePhase', signature(object='data.frame'),
  function(object,xlim=c(0,ceiling(2*max(object$stock,  na.rm=TRUE))/2),
                  ylim=c(0,ceiling(2*max(object$harvest,na.rm=T))/2)){
    
       invisible(kobePhaseFn(object,xlim,ylim))})

