#' @title Plot an estimated CUB model
#' @description Plotting facility for the CUB estimation of ordinal responses. 
#' @aliases cubvisual
#' @usage cubvisual(m, ordinal, caption, labelpoint="estim", xlim, ylim)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param caption Characters string for the plot caption (default is "CUB models parameter space")
#' @param labelpoint Characters string describing the label to give to the point estimates (default is "estim")
#' @param xlim Numeric vectors of length 2, giving the x coordinate range (if missing, xlim=c(0,1))
#' @param ylim Numeric vectors of length 2, giving the y coordinate range (if missing, ylim=c(0,1))
#' @details The estimation is performed via \code{\link{cubforsim}}. It represents an estimated CUB model as a point
#'  in the parameter space with some useful options. If necessary, other estimated models may be added with
#' the standard commands of the R environment (as points(.), for instance).
#' @return A plot of the estimated parameter vector \eqn{(\pi, \xi)} as a point in the parameter space
#' @keywords device
#' @export cubvisual
#' @import graphics
#' @examples
#' data(univer)
#' m<-7
#' ordinal<-univer[,12] 
#' cubvisual(m, ordinal)

cubvisual <-
function(m,ordinal,caption,labelpoint="estim",xlim,ylim){
  
  if (missing(xlim)){
    xlim<-c(0,1)
  }
  if (missing(ylim)){
    ylim<-c(0,1)
  }
  if (missing(caption)){
    caption<-"CUB models parameter space"
  }
  
  stimacub<-cubforsim(m, ordinal, maxiter = 500, toler = 1e-06)
    #cub00(m,ordinal,maxiter=500,toler=1e-6,makeplot=FALSE);
  param<-stimacub$estimates; pai<-param[1];csi<-param[2];
  plot(1-pai,1-csi,main=caption,las=1,pch=19,cex=1.2,xlim=xlim,ylim=ylim,
       col="blue",
       xlab=expression(paste("Uncertainty  ", (1-pi))),
       ylab=expression(paste("Feeling  ", (1-xi))));
  text(1-pai,1-csi,labels=labelpoint,font=4,pos=1,offset=0.4,cex=0.8)
}
