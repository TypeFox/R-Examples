plot.Simpson <-
function(x,...)
{
  alldata <- x$alldata
  Allbeta <- x$Allbeta
  Allint <- x$Allint
  Nclusters <- x$Nclusters
  namex <- x$namex
  namey <- x$namey
  
  #plotting all clusters and drawing regression lines

  plot(alldata[,1],alldata[,2], col = alldata[,3]+1,pch=alldata[,3],xlab=namex,ylab=namey)
  Allbeta=as.numeric(Allbeta)
  Allint=as.numeric(Allint)
  groupreg=lm(alldata[,2]~alldata[,1])
  groupint=groupreg$coefficients[1]
  groupbeta=groupreg$coefficients[2]
  for(j in 1:Nclusters)
  {
    a=Allint[j]
    b=Allbeta[j]
    abline(a,b,lty=j,col=j+1)
  }
  abline(groupint,groupbeta,lwd=3,col=1)
  xpd <- par("xpd")
  par(xpd=TRUE)
  legend(mean(par("usr")[c(1,2)]),par("usr")[4],paste("cluster", 1:Nclusters),col= 1:Nclusters+1, pch=1:Nclusters,bty="n",horiz=TRUE,xjust=0.5,yjust=0,cex=1)
  par(xpd=xpd)
  }
