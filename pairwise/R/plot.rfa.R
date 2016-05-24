#' @export plot.rfa
#' @title S3 Plotting Rasch Residual Factor Analysis
#' @description S3 plotting Method for object of class\code{"rfa"}
#' @param x object of class\code{"rfa"}
#' @param com an integer giving the number of the principal component used for plotting
#' @param ra either the character \code{"auto"} (default) or an numeric, defining the (logit) range for x-axis
#' @param main see \code{\link{plot}}
#' @param labels a character vector specifying the plotting pattern to use. see \code{\link{text}}. At default the itemnames are used.
#' @param xlab see \code{\link{plot}}
#' @param ylab see \code{\link{plot}}
#' @param srt see \code{\link{text}} or \code{\link{par}}
#' @param cex.axis see \code{\link{plot}}
#' @param cex.text see argument \code{cex} in function \code{\link{text}}
#' @param col.text see argument \code{col} in function \code{\link{text}}
#' @param ... other parameters passed through.

########################### hier die plot method für rfa #############################
plot.rfa<-function(x, com=1, ra="auto", main=NULL , labels=NULL, xlab="logits", ylab="loadings", srt=0, cex.axis = 0.8, cex.text = 0.8, col.text = NULL, ...){
  if(length(main)==0){main<-paste(com, " Component for ","\n", deparse(substitute(x)),sep="")}
  if(length(labels)==0){labels <- rownames(x$pca$loadings)}
  #names(x$pers_obj$pair$sigma)
  
  bereich <- ra
  
  yy <- x$pca$loadings[,com]
  if(x$transposed==FALSE){xx <- x$pers_obj$pair$sigma} # added: 17-12-2014 
  if(x$transposed==TRUE){xx <- x$pers_obj$pers$WLE} # added: 17-12-2014
  
   ##### plotingrange festlegen mit leerplot
  ## automatische x achsen skalierung
  if((bereich)[1]=="auto"){
    xx1<-floor(min(xx))
    xx2<-ceiling(max(xx))
  }
  ## feste vorgegebene x achsen skalierung
  if(class(bereich)=="numeric"){
    xx1<- -bereich
    xx2<- bereich
  }
  
  yauto <- c(((floor(min(yy)*10))/10), ((ceiling(max(yy)*10))/10))
  
  # empty ploting range
  plot(x=c(xx1,xx2), y=yauto , type="n",xlab=xlab,ylab=ylab,bty="n",main=main, cex.axis=cex.axis)
  text(x=xx,y=yy,labels=labels, cex=cex.text, srt=srt, col=col.text, ...)
}
  