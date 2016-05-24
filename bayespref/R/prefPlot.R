prefPlot <-
function(prefres=NULL,burn=0,ind=TRUE,pop=TRUE,dist=FALSE,ymax=5,
                   xmin=0,xmax=1,dadj=2,colors=NULL,leg=FALSE,lx=0.8,ly=4,
                   catname=NULL,ps=FALSE,file="./prefplot.eps"){

  allcolors<-brewer.pal(n=8,"Accent") # default colors
  mcmc<-dim(prefres$PopPref)[2]
  nCat<-dim(prefres$PopPref)[1]
  nInd<-dim(prefres$IndPref)[1]
  if (is.null(colors)==TRUE){ # user didn't supply colors use mine
    colors<-rep(allcolors,4)[1:nCat] #reuse colors if more than 8 cats, this is good for up to 32
  }
  if (ps==TRUE){
    postscript(file=file,width=6,height=6)
  }
  
  plot(0:1,0:1,type='n',ylim=c(0,ymax),xlim=c(xmin,xmax),xlab="Preference",
       ylab="Probability",cex.lab=1.2,cex.axis=1.1)

  if (ind==TRUE){
    for (i in 1:nInd){
      for (j in 1:nCat){
        lines(density(prefres$IndPref[i,j,(burn+1):mcmc]),col=colors[j],lty=2)
      }
    }
  }

  if (pop==TRUE){
    for (j in 1:nCat){
      lines(density(prefres$PopPref[j,(burn+1):mcmc],adj=dadj),col=colors[j],lty=1,lwd=3)#JF add lwd=3
    }
  }

  if (dist==TRUE){
    alphas<-numeric(nCat)
    for (j in 1:nCat){
      alphas[j]<-median(prefres$PopPref[j,(burn+1):mcmc]) * median(prefres$PopVar[(burn+1):mcmc])
    }
    x<-rdirichlet(10000,alphas)
    for (j in 1:nCat){
      lines(density(x[,j]),col=colors[j],lty=1,lwd=4)
    }
  }

  if (leg==TRUE){
    catname2<-numeric(nCat * 2)
    for (i in 1:nCat){
      catname2[i]<-paste(catname[i],"pop")
    }
    for (i in (nCat+1):(2*nCat)){
      catname2[i]<-paste(catname[(i-nCat)],"ind")
    }
    if(ind==TRUE){legend(lx,ly,catname2,col=c(colors[1:nCat],colors[1:nCat]),lty=c(rep(1,nCat),rep(2,nCat)),
           lwd=c(rep(2,nCat),rep(1,nCat)))}
    if(ind==FALSE){legend(lx,ly,catname2[1:nCat],col=c(colors[1:nCat],colors[1:nCat]),lty=c(rep(1,nCat),rep(2,nCat)),
           lwd=c(rep(2,nCat)))}
       
  }
  if (ps==TRUE){
    dev.off()
  }
}

