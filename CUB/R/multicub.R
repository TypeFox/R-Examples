#' @title Joint plot of estimated CUB models in the parameter space
#' @description Return a plot of estimated CUB models represented as points in the parameter space.
#' @aliases multicub
#' @usage multicub(matord, m, labelpoints=as.character(1:ncol(matord)), 
#' caption="CUB models", colours="black", symbols=19, 
#' thickness=1.5, xwidth=c(0,1), ywidth=c(0,1))
#' @export multicub
#' @param matord Matrix of ordinal data arranged per columns
#' @param m Number of ordinal categories
#' @param labelpoints Character strings indicating the labels for the estimated models as points
#'  in the parameter space (default is the corresponding index of column in matord)
#' @param caption Character string indicating the plot title (default is "CUB models")
#' @param colours Indicate the colours for the plotted points
#' @param symbols Indicate the symbols used to represent estimated CUB models as points in the
#'  parameter space
#' @param thickness Indicate the thickness
#' @param xwidth Indicate the width of the abscissa axis
#' @param ywidth Indicate the width of the ordinate axis
#' @keywords device
#' @examples
#' data(univer)
#' matord<-univer[,8:12]
#' m<-7
#' multicub(matord, m, labelpoints=as.character(1:ncol(matord)), 
#'         caption ="CUB models", colours="black", symbols=19,
#'         thickness=1.5, xwidth=c(0,1), ywidth=c(0,1))

multicub <-
function(matord,m,labelpoints=as.character(1:ncol(matord)),
                   caption="CUB models",colours="black", symbols=19,
                   thickness=1.5,xwidth=c(0,1),ywidth=c(0,1)){
  k<-ncol(matord)
  vettpai<-vettcsi<-rep(NA,k);
  for(j in 1:k){
    stimacub <- cubforsim(m,matord[,j])
    #stimacub <- cub00(m,matord[,j],maxiter=500,toler=1e-6,makeplot=FALSE);
    param<-stimacub$estimates
    vettpai[j]<-param[1]; vettcsi[j]<-param[2];
  }
  plot(1-vettpai,1-vettcsi,main=caption,cex=1.2,cex.main=1, font.lab=4,cex.lab=1,
       pch=symbols, col=colours,lwd=thickness,las=1,
       xlim=xwidth,ylim=ywidth,
       xlab=expression(paste("Uncertainty  ", (1-pi))),
       ylab=expression(paste("Feeling  ", (1-xi))));
  text(1-vettpai,1-vettcsi,labels=labelpoints,pos=3,offset=0.5,font=4,cex=0.9)
}
