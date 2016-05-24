`ColorDendrogram` <-
function(hc, y, main="", branchlength=.7, labels=NULL, xlab=NULL, sub=NULL,ylab="",cex.main=NULL){
  if(is.null(labels)) labels <- rep("", length(y))
  plot(hc,hang=.2, main=main, labels=labels, xlab=xlab, sub=sub,ylab=ylab,cex.main=cex.main) #plclust
  y <- y[hc$ord]
  if(is.numeric(y)){
    y <- y+1
    y[y==7] <- "orange"
  } 
  for(i in 1:length(hc$ord)){
    o=hc$merge[,1]==  -hc$ord[i] | hc$merge[,2]== -hc$ord[i]
    segments(i,hc$he[o]-branchlength,i,hc$he[o],col=y[i])
  }
}

