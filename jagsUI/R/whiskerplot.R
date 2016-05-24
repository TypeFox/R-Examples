
whiskerplot <- function(x,parameters,quantiles=c(0.025,0.975),zeroline=TRUE){
  if(class(x)!="jagsUI"){stop('Requires jagsUI object as input')}
  devAskNewPage(ask=FALSE)
  
  #Generate a list of all specified output parameters
  #Expand from shorthand if necessary
  parameters <- translate.params(x,parameters)
  
  n <- length(parameters)
  
  xstructure <- c(1:n)
  
  qs <- function(x,y){as.numeric(quantile(x,y))}
  
  means <- tops <- bottoms <-ymin <- ymax <- vector(length=n)
  for (i in 1:n){
    hold <- unlist(x$samples[,parameters[i]])
    means[i] <- mean(hold)
    tops[i] <- qs(hold,quantiles[2])
    bottoms[i] <- qs(hold,quantiles[1])   
  }
  
  ymin <- min(bottoms)
  ymax <- max(tops)
  
  plot(xstructure,means,xaxt="n",ylim=c(ymin,ymax),xlim=c(0,n+1),xlab="Parameters",ylab="Parameter Values",pch=19,cex=1.5,
       main=paste('Whisker plot, quantiles (',quantiles[1],' - ',quantiles[2],')',sep=""))
  axis(side=1, at=c(1:n), labels=parameters)
  box()
  
  if(zeroline){abline(h=0)}
  
  for (i in 1:n){
    segments(x0=xstructure[i],y0=bottoms[i],x1=xstructure[i],y1=tops[i], lwd=2)
    segments(x0=xstructure[i]-0.2,y0=bottoms[i],x1=xstructure[i]+0.2,y1=bottoms[i])
    segments(x0=xstructure[i]-0.2,y0=tops[i],x1=xstructure[i]+0.2,y1=tops[i])
  }
 
}