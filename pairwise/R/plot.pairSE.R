#' @export plot.pairSE
#' @title S3 Plotting Thustonian Thresholds with SE
#' @description S3 plotting method for object of class\code{"pairSE"}
#' @param x object of class\code{"pairSE"}
#' @param sortdif logical wether to order items by difficulty
#' @param ra either the character \code{"auto"} (default) or an numeric, defining the (logit) range for y-axis
#' @param ci numeric defining confidence intervall for point estimator
#' @param main see \code{\link{plot}}
#' @param col.lines vector of colors for threshold profile lines
#' @param col.error vector of colors for error bars
#' @param type see \code{\link{plot}}
#' @param xlab see \code{\link{plot}}
#' @param ylab see \code{\link{plot}}
#' @param pch see \code{\link{plot}}
#' @param las see \code{\link{plot}}
#' @param cex.axis see \code{\link{plot}}
#' @param ... other parameters passed to plot
########################### hier die plot method für pairSE #############################
#ci=2; sortdif=FALSE; ra="auto"; main=NULL; col.lines=1:(dim(x$parameter)[2]-1); col.error=1:(dim(x$parameter)[2]-1); type="b";xlab="items"; ylab="logits"; pch=20; las=3; cex.axis = 0.8

#plot.pairSE<-function(x, ci=2, sortdif=FALSE, ra="auto", main=NULL, col.lines=1:(dim(x$parameter)[2]-1), col.error=1:(dim(x$parameter)[2]-1), type="b",xlab="items", ylab="logits", pch=20, las=3, cex.axis = 0.8, ...){
plot.pairSE<-function(x, ci=2, sortdif=FALSE, ra="auto", main=NULL, col.lines=1:(dim(x$threshold)[2]), col.error=1:(dim(x$threshold)[2]), type="b",xlab="items", ylab="logits", pch=20, las=3, cex.axis = 0.8, ...){
  
  if(length(main)==0){main<-deparse(substitute(x))}
  
  bereich <- ra
  
  if(sortdif==TRUE){
    #sorter <- order(x$parameter[,"sigma"])
    sorter <- order(x$sigma)
    # x$parameter <- x$parameter[sorter,]
    x$SE <- x$SE[sorter,]
    x$SEsigma <- x$SEsigma[sorter]
    x$threshold <- x$threshold[sorter,] 
    x$sigma <- x$sigma[sorter] 
    cat("(ordered by location) \n")
  }
  #ende der sortierung
  # thresholds<-as.matrix(x$parameter[,1:(dim(x$parameter)[2]-1)])
  thresholds <- x$threshold
  # SEthresholds<-as.matrix(x$SE[,1:(dim(x$SE)[2]-1)])
  SEthresholds <- x$SE
  if(dim(thresholds)[2]==1){colnames(thresholds)="sigma"}
  thresholds->logit
  maxLen <- dim(logit)[2]; nitem <- dim(logit)[1]

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
  plot(c((1-.5),(nitem +.5) ),c(y1,y2), type="n",xaxt="n",xlab=xlab,ylab=ylab,bty="n",main=main)
  # plotting lines
  matplot(logit,add=TRUE,pch=pch,type=type, xlab=xlab, xaxt="n", main=main, col = col.lines)#, ...
  axis(1, 1:(nitem), labels=c(rownames(logit)),las=las, cex.axis=cex.axis)#, ... 
  mtext(text=colnames(logit), side = 4, at = colMeans(logit,na.rm=TRUE), padj = NA, cex = cex.axis, col = col.lines,las=1)
  
  #for (i in 1:(dim(x$parameter)[2]-1)){
  for (i in 1:(dim(x$threshold)[2])){
    #segments( 1:(dim(x$parameter)[1]), thresholds[,i]+SEthresholds[,i]*ci, 1:(dim(x$parameter)[1]), thresholds[,i]-SEthresholds[,i]*ci ,col=col.error[i])#,... 
    segments( 1:(dim(x$threshold)[1]), thresholds[,i]+SEthresholds[,i]*ci, 1:(dim(x$threshold)[1]), thresholds[,i]-SEthresholds[,i]*ci ,col=col.error[i])#,... 
    
    #segments( 1:(dim(x$parameter)[1])-((.15*(SEthresholds[,i]!=0))*ci), thresholds[,i]+SEthresholds[,i]*ci, 1:(dim(x$parameter)[1])+((.15*(SEthresholds[,i]!=0))*ci), thresholds[,i]+SEthresholds[,i]*ci ,col=col.error[i])#,...
    segments( 1:(dim(x$threshold)[1])-((.15*(SEthresholds[,i]!=0))*ci), thresholds[,i]+SEthresholds[,i]*ci, 1:(dim(x$threshold)[1])+((.15*(SEthresholds[,i]!=0))*ci), thresholds[,i]+SEthresholds[,i]*ci ,col=col.error[i])#,... 
    
    #segments( 1:(dim(x$parameter)[1])-((.15*(SEthresholds[,i]!=0))*ci), thresholds[,i]-SEthresholds[,i]*ci, 1:(dim(x$parameter)[1])+((.15*(SEthresholds[,i]!=0))*ci), thresholds[,i]-SEthresholds[,i]*ci ,col=col.error[i])#,...  
    segments( 1:(dim(x$threshold)[1])-((.15*(SEthresholds[,i]!=0))*ci), thresholds[,i]-SEthresholds[,i]*ci, 1:(dim(x$threshold)[1])+((.15*(SEthresholds[,i]!=0))*ci), thresholds[,i]-SEthresholds[,i]*ci ,col=col.error[i])#,... 
  }
  par(op) # reset graphics setting
  
}
