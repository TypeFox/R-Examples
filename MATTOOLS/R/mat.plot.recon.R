 
## The function is currently defined as
mat.plot.recon<-function (inrec,colAgDp=7,inCritVal) 
{

tmethod=attributes(inrec)$weightmethod
reccol = which(attributes(inrec)$names == "recname")

if(tmethod != "none")
        {

        windows(height=11,width=8.5)
        par(mfrow=c(2,1))
        xnt=c(inrec[,c(reccol)],inrec[,c(reccol+1)],inrec[,c(reccol+2)])

        plot(c(min(inrec[,colAgDp],na.rm=T),max(inrec[,colAgDp],na.rm=T)),c(min(xnt,na.rm=T),max(xnt,na.rm=T)),type="n",xlab="Age/Depth",ylab="Environmental Value")

        lines(inrec[,colAgDp],inrec[,reccol])

        lines(inrec[,colAgDp],inrec[,reccol+1],col="blue")

        lines(inrec[,colAgDp],inrec[,reccol+2],col="blue")

        }

if(!missing(inCritVal)){

                plot(inrec[,colAgDp],inrec[,reccol+3],type="l",xlab="Age/Depth",ylab="Dissimilarity")

                abline(h=inCritVal$cutoffs$y,lwd=2,lty=3,col="red")

                textloc = max(inrec[,colAgDp],na.rm=T)-((max(inrec[,colAgDp],na.rm=T)-min(inrec[,colAgDp],na.rm=T))/50)

                text(textloc,inCritVal$cutoffs$y,round(inCritVal$cutoffs$x,3),pos=3,col="red")

        }

  }



