#---------------------------------------------------------------#
# restScore  Author Niels G Waller
# Latest revision January 27, 2016
#---------------------------------------------------------------#
restScore<- function(data,  item, NCuts=10){

	Nitems<-ncol(data)
   tot.score<-apply(data,1,sum)
   pmat<-rep(0,NCuts)
   rest.score<-tot.score - data[,item]
 
    probs<-unlist(lapply(split(data[,item],cut(rest.score,NCuts)),mean))
    group.N<-unlist(lapply(split(data[,item],cut(rest.score,NCuts)),length))

    se.p<- sqrt(( probs * (1-probs))/group.N) 
    plot(1:NCuts,probs,axes=FALSE,ylim=c(0,1),
        xlab="Rest Scores",
        ylab="Probability  of  a  Keyed  Response",
        pch=16,col="blue",cex=1.5,
        lwd=3,
        main=paste("Item ",item, sep="")) 

    # bound limits of std error bars to (0,1)
    upperProbs <- probs + 1.96*se.p
    upperProbs[upperProbs>1]<-1
    lowerProbs <- probs - 1.96*se.p
    lowerProbs[lowerProbs < 0] <- 0

    for (i in 1:NCuts) { 
        lines( c((1:NCuts)[i],(1:NCuts)[i]),c(probs[i],upperProbs[i]),lwd=3, lty=4) 
        }
    for (i in 1:NCuts) { 
        lines( c((1:NCuts)[i],(1:NCuts)[i]),c(probs[i],lowerProbs[i]),lwd=3, lty=4) 
        }         

    
    axis(side=1,at=1:NCuts, pos=0,
      labels=names(table(cut(rest.score,NCuts))),
      cex.axis=.68,
      lwd=3)
    
    
    axis(side=2,lwd=3,pos= .75,cex.axis=.68)
    abline(h=1,lwd=3)
    abline(h=0,lwd=3)

    points(1:NCuts,probs,type="l",lwd=3, pch=1, 
      cex=1.5,
      ylim=c(0,1))
   
   # draw lines at upper and lower probability bounds
    if((probs[1]-1.96*se.p[1])>0)
     abline(h=probs[1]-1.96*se.p[1],col="lightgrey",lty=3,lwd=3)
    if((probs[NCuts]+1.96*se.p[NCuts])<1)
     abline(h=probs[NCuts]+1.96*se.p[NCuts],col="grey",lty=3,lwd=3)
  
    invisible(list(item = item, bins = table(cut(rest.score,NCuts)),binProb = probs ))
}    


