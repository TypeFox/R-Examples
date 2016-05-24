#' @export plot.grm
#' @title S3 Plotting Graphical Model Check
#' @description S3 plotting Method for object of class\code{"grm"}
#' @param x object of class\code{"grm"}
#' @param xymin optional lower limit for xy-axis
#' @param xymax optional upper limit for xy-axis
#' 
#' @param ci numeric defining confidence intervall for point estimator
#' @param main see \code{\link{plot}}
#' @param col.error vector of colors for error bars
#' @param col.diag color for the diagonal of the plot
#' @param itemNames logical wether to plot itemnames
#' @param cex.names magnification factor for itemnames
#' @param type see \code{\link{plot}}
#' @param xlab see \code{\link{plot}}
#' @param ylab see \code{\link{plot}}
#' @param pch see \code{\link{plot}}
#' @param las see \code{\link{plot}}
#' @param cex.axis see \code{\link{plot}}
#' @param ... other parameters passed to plot
########################### hier die plot method für grm #############################
plot.grm<-function(x, xymin=NULL, xymax=NULL, ci=2, main=NULL, col.error="blue", col.diag="red", itemNames=TRUE, cex.names=.8, type="b", xlab=NULL, ylab=NULL, pch=43, las=3, cex.axis = 0.5, ...){  
  
  if(length(main)==0){main<-deparse(substitute(x))}
  
  #sonderfall 2 subsamples
  if (length(x)==2){
    subsamp_names <- names(x)
    
#    Itemnames<-names(x[[1]]$threshold[,dim(x[[2]]$threshold)[2]])
    Itemnames <- rownames(x[[1]]$threshold)
#     X<-x[[subsamp_names[1]]]$threshold[,dim(x[[ subsamp_names[1] ]]$threshold)[2]]
#     Y<-x[[subsamp_names[2]]]$parameter[,dim(x[[ subsamp_names[2] ]]$parameter)[2]]
    X<-x[[subsamp_names[1]]]$sigma
    Y<-x[[subsamp_names[2]]]$sigma
    
#    XS<-x[[ subsamp_names[1] ]]$SE[,dim(x[[ subsamp_names[1] ]]$SE)[2]]
#    YS<-x[[ subsamp_names[2] ]]$SE[,dim(x[[ subsamp_names[2] ]]$SE)[2]]
    XS<-x[[ subsamp_names[1] ]]$SEsigma
    YS<-x[[ subsamp_names[2] ]]$SEsigma

    ##### plotingrange festlegen mit leerplot
    if(length(xymax)==0){xymax<-(round((max(c(max(X),max(Y))) + 3*(max(c(max(XS),max(YS)))))*10))/10}
    if(length(xymin)==0){xymin<-(round((min(c(min(X),min(Y))) - 3*(max(c(max(XS),max(YS)))))*10))/10}
    
    xx<-c(xymin,xymax); yy<-c(xymin,xymax)
    if (length(xlab)==0) {xlab <- subsamp_names[1]}
    if (length(ylab)==0) {ylab <- subsamp_names[2]}
    
    plot(x=xx,y=yy,type="n",bty="n", main=main, xlab=xlab, ylab=ylab, ...)
    # hilfsfunktion elipse
    eli<-function(x.cent,y.cent,xb,yh, ...){
      # plotten einer elipse
      nseg=360
      xx <- x.cent + xb*cos( seq(0,2*pi, length.out=nseg) )
      yy <- y.cent + yh*sin( seq(0,2*pi, length.out=nseg) )
      lines(xx,yy,col=col.error,...) 
    }
    ##############
    ##### plotten der grafik    
    if (itemNames==TRUE){pch=""}
    if (itemNames==TRUE){text(X,Y,Itemnames,cex=cex.names,...)}
    points(X,Y,pch=pch,...) 
    abline(0,1,col=col.diag,...) 
    for (i in 1: length(X)){ eli(X[i],Y[i],XS[i]*ci,YS[i]*ci) }
}
  
  if (length(x)!=2){cat("actualy no plotting method for", length(x) ,"subsample(s) available","\n")  }
  
}
  