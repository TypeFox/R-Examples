#' @export plot.pair
#' @title S3 Plotting Thurstonian Thresholds
#' @description S3 plotting Method for object of class\code{"pair"}
#' @param x object of class\code{"pair"}
#' @param sortdif logical wether to order items by difficulty
#' @param ra either the character \code{"auto"} (default) or an numeric, defining the (logit) range for y-axis
#' @param main see \code{\link{plot}}
#' @param col.lines vector of colors for threshold profile lines
#' @param type see \code{\link{plot}}
#' @param xlab see \code{\link{plot}}
#' @param ylab see \code{\link{plot}}
#' @param pch see \code{\link{plot}}
#' @param las see \code{\link{plot}}
#' @param cex.axis see \code{\link{plot}}

#' @param ... other parameters passed to plot

########################### hier die plot method für pair #############################
plot.pair<-function(x, sortdif=FALSE, ra="auto", main=NULL , col.lines=(1:dim(x$threshold)[2]), type="b", xlab="items", ylab="logits", pch=(1:dim(x$threshold)[2]) , las=3, cex.axis = 0.8, ...){
  if(length(main)==0){main<-deparse(substitute(x))}
  
  bereich <- ra
  
  if(sortdif==TRUE){
    threshold <- x$threshold
    #sb <- x$sb
    sigma <- x$sigma
    #####
    threshold <- threshold[order(sigma), ]
    #sb <- sb[order(sigma), ]
    sigma <- sort(sigma)
    x<-list(threshold=threshold,sigma=sigma)
    class(x)<-c("pair", "list")
    cat("(ordered by location) \n")
  }
  
  x$threshold->logit
  maxLen <- dim(logit)[2]
  colnames(logit)<-paste("threshold",1:maxLen,sep=".")
  
  op <- par(mar = c(5,4,4,6),bty="n",oma=c(0, 0, 0, 0) ) # set graphics options
  
  ##### plotingrange festlegen mit leerplot
  ## automatische y achsen skalierung
  if((bereich)[1]=="auto"){
    y1<-floor(min(logit,na.rm=TRUE))
    y2<-ceiling(max(logit,na.rm=TRUE))
  }
  ## feste vorgegebene y achsen skalierung
  if(class(bereich)=="numeric"){
    y1<- -bereich
    y2<- bereich
  }
  # empty ploting range
  plot(c(1,(dim(logit)[1])),c(y1,y2), type="n",xaxt="n",xlab=xlab,ylab=ylab,bty="n",main=main)
  matplot(logit,add=TRUE,type=type, xlab=xlab, xaxt="n", main=main, col = col.lines)#, ...
  axis(1, 1:(dim(logit)[1]), labels=c(rownames(logit)),las=las, cex.axis=cex.axis)#, ...
  axis(4, at=colMeans(logit,na.rm=T) , labels=colnames(logit) , outer = F, las=1,cex.axis=.8 ,lwd=0  )

#text(x=rep(dim(logit)[1],dim(logit)[2]), y = seq((y1*.9),(y2*.9), length.out=dim(logit)[2]), labels=colnames(logit),pos=4,col=col.lines,cex=cex.axis)

par(op) # reset graphics setting
}
  