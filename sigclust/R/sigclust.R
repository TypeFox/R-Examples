# function "cluster" perform 2-means clustering and then calculate clusting
# index

#input data matrix x of size n by p: n is sample size (no. of subjects) and p
#is dimension (no. of genes)

.cluster<-function(x, n, p){
  if(n>1){
    x<-as.matrix(x)
                                        #check the dimension of x to match
                                        #n and p
    if(dim(x)[1]==n & dim(x)[2]==p){	#run 2-means on x
      clust<-kmeans(x,2)
      withinsum<-sum(clust$withinss)
      meanp<-colMeans(x)
      tx<-t(x)
      txdiff<-tx-meanp
      totalsum<-sum(txdiff^2)
      cindex<-withinsum/totalsum
      list(clust=clust, cindex=cindex)
    }else{
      print("Wrong size of matrix x!")
      return(0)
    }
  }else{
    print("Only one sample left, no need for clustering!")
    return(0)
  }
  
}

.sigclustcovest <- function(vsampeigv,sig2b){
  d <- length(vsampeigv)
                                        #Check have some eigenvalues < sig2b
  vtaucand <- vsampeigv - sig2b
                                        #find threshold to preserve power
  which <- which(vtaucand<=0)
  icut <- which[1] - 1
  powertail <- sum(vsampeigv[(icut+1):d])
  power2shift <- sig2b*(d-icut) - powertail
  vi <- c(1:icut)
  vcumtaucand <- sort(cumsum(sort(vtaucand[vi])),decreasing=TRUE)
  vpowershifted <- (vi-1)*vtaucand[vi] + vcumtaucand
  flag <- vpowershifted < power2shift
  if(sum(flag)==0){
    itau <- 0;
  }else{
    which <- which(flag>0)
    itau <- which[1]
  }
  if(itau==1){
    powerprop <- power2shift/vpowershifted
    tau <- powerprop*vtaucand[1]
  }else if(itau==0){
    powerprop <- power2shift/vpowershifted[icut] 
    tau <- powerprop*vtaucand[icut] 
  }else{
    powerprop <- (power2shift-vpowershifted[itau])/
      (vpowershifted[itau-1]-vpowershifted[itau]) 
    tau <- vtaucand[itau] + powerprop*(vtaucand[itau-1] - vtaucand[itau]) 
  }
  veigvest <- vsampeigv - tau 
  flag <- veigvest > sig2b 
  veigvest <- flag*veigvest + (1-flag)*(sig2b*rep(1,d))
  list(veigvest=veigvest,tau=tau)
}

  
                                        #calculate variances of normal for
                                        #simulation: eigen values of pca and
                                        #background variance.

.vareigen<-function(x,n,p,icovest){	
                                        #check the dimension of x to
                                        #match n and p
  if(dim(x)[1]==n & dim(x)[2]==p){
    mad1<-mad(x)
    simbackvar<-mad1^2
                                        # Jan. 23, 07; replace eigen by
                                        #svd to save memory
    xcov<-cov(x)
    xeig<-eigen(xcov, symmetric=TRUE, only.values =TRUE)
    veigval<-xeig$values
    vsimeigval<-xeig$values
    
                                        #      avgx<-t(t(x)-colMeans(x))
                                        #	dv<-svd(avgx)$d
                                        #      veigval<-dv^2/(n-1)
                                        #	vsimeigval<-veigval

    if(icovest==1){
      taub = 0
      tauu <- .sigclustcovest(veigval,simbackvar)$tau
      etau = (tauu-taub)/100
      ids = rep(0,100)
      for(i in 1:100){
        taus = taub + (i-1)*etau
        eigval.temp <- veigval - taus
        eigval.temp[eigval.temp<simbackvar] <- simbackvar
        ids[i] <- eigval.temp[1]/sum(eigval.temp)
      }
      tau <- taub + (which.max(ids)-1)*etau
      vsimeigval <- veigval - tau
      vsimeigval[vsimeigval<simbackvar] <- simbackvar
      #vsimeigval <- .sigclustcovest(veigval,simbackvar)$veigvest
    }

    if(icovest==2){
      vsimeigval[veigval<0] <- 0
    }
    
    if(icovest==3){  # Use original background noise thresholded estimate
                                        # (from Liu, et al, JASA paper)
      vsimeigval[veigval<simbackvar] <- simbackvar
    }
    list(veigval=veigval, simbackvar=simbackvar, vsimeigval=vsimeigval)
  }else{
    print("Wrong size of matrix x!")
    return(0)
  }
}


#give null variances (eigenvalues), simulate normal vectors and compute
#cluster index

.simnull<-function(vsimeigval, n, p){
  simnorm<-matrix(0, n, p)
  for(i in 1:n){
    simnorm[i,]<-rnorm(p, sd=sqrt(vsimeigval))
  }
  simclust<-.cluster(simnorm, n, p)
  list(cindex=simclust$cindex)
}

#repeat simulation "nsim" times and obtain distribution of clustering index;
#then calculate p-value

#labflag is an indicator variable specifying if the p-value is for assigned
#cluster or use 2-means; for assiged cluster labflag=1, otherwise labflag=0 if
#labflag=1, label is a vector giving the assigment of cluster (1,2) if
#labflag=0, label is 0

setClass("sigclust",
         representation(raw.data="matrix",
                        veigval="vector",
                        vsimeigval="vector",
                        simbackvar="numeric",
                        icovest="numeric",
                        nsim="numeric",
                        simcindex="vector",
                        pval="numeric",
                        pvalnorm="numeric",
                        xcindex="numeric"))

sigclust <- function(x,nsim,nrep=1,labflag=0,label=0,icovest=1){	
                                        #check the dimension of x to
                                        #match n and p
  n <- dim(x)[1]
  p <- dim(x)[2]
  if(n>1){
    x<-as.matrix(x)
    if(labflag==0){
      xclust<-.cluster(x, n, p)
                                        #clust.0<-cluster(x,n,p)
      for(i in 1:nrep){
        clust.temp<-.cluster(x,n,p)
        if(clust.temp$cindex<xclust$cindex)
          xclust<-clust.temp
        xcindex<-xclust$cindex
      }
    }
    xvareigen<-.vareigen(x,n,p,icovest)
    simcindex<-rep(0,nsim)
    for(i in 1:nsim){
      xsim<-.simnull(xvareigen$vsimeigval, n, p)
      simcindex[i]<-xsim$cindex
    }
    if(labflag==0){
      index<-(simcindex<=xclust$cindex)
      mindex<-mean(simcindex)
      sindex<-sd(simcindex)
      pval<-sum(index)/nsim
      pvalnorm<-pnorm(xclust$cindex,mindex,sindex)
    }
    if(labflag==1){
                                        #first calcuate the cluster index for x
      meanp1<-colMeans(x[label==1,])
      txdiff1<-t(x[label==1,])-meanp1
      meanp2<-colMeans(x[label==2,])
      txdiff2<-t(x[label==2,])-meanp2
      withinsum<-sum(txdiff1^2)+sum(txdiff2^2)
      meanp<-colMeans(x)
      tx<-t(x)
      txdiff<-tx-meanp
      totalsum<-sum(txdiff^2)
      cindexlab<-withinsum/totalsum
      index<-(simcindex<=cindexlab)
      mindex<-mean(simcindex)
      sindex<-sd(simcindex)
      pval<-sum(index)/nsim
      pvalnorm<-pnorm(cindexlab,mindex,sindex)
      xcindex<-cindexlab
      
    }
    
    return(new("sigclust",raw.data=x,
               veigval=xvareigen$veigval,
               vsimeigval=xvareigen$vsimeigval,
               simbackvar=xvareigen$simbackvar,
               icovest=icovest,
               nsim=nsim,
               simcindex=simcindex,
               pval=pval,
               pvalnorm=pvalnorm,
               xcindex=xcindex))
  }else{
    print("Only one sample left, no need for clustering!")
    return(0)
  }
  
}

#plot p-value with cluster index and the density of simulated cluster indexes.

.plot.sigclust <- function(sigclust,arg="all",...){
  
  raw.data <- sigclust@raw.data
  veigval <- sigclust@veigval
  simbackvar <- sigclust@simbackvar
  vsimeigval <- sigclust@vsimeigval
  icovest <- sigclust@icovest
  nsim <- sigclust@nsim
  simcindex <- sigclust@simcindex
  pval <- sigclust@pval
  pvalnorm <- sigclust@pvalnorm
  xcindex <- sigclust@xcindex

  n <- dim(raw.data)[1]
  d <- dim(raw.data)[2]

  #Background Standard Deviation Diagnostic Plot

  overlay.x <- as.vector(raw.data)
  mean <- mean(overlay.x)
  sd <- sd(overlay.x)
  median <- median(overlay.x)
  mad <- mad(overlay.x)
  ntot <- length(overlay.x)
  maxnol <- 5000

  if(arg=="background"|arg=="all"){
    par(mfrow=c(1,1))
    par(mar=c(5,4,4,2)+0.1)
    nused <- maxnol
    denraw <- density(overlay.x)
    if(ntot>maxnol){
      overlay.x <- overlay.x[sample(c(1:ntot),maxnol)]
    }else{
      nused <- ntot
    }
    
    xmin <- min(denraw$x)
    xmax <- max(denraw$x)
    ymin <- min(denraw$y)
    ymax <- max(denraw$y)
    
    overlay.y <- ymin + (0.15+0.5*runif(nused))*(ymax - ymin)
    plot(denraw,xlim=range(overlay.x),ylim=range(denraw$y),
         col="blue",xlab="",main="",lwd=3,...)
    xgrid <- seq(xmin,xmax,by=0.0025*(xmax-xmin))
    normden <- dnorm(xgrid,mean=median,sd=sd)
    lines(xgrid,normden,col="red",lwd=3)
    points(overlay.x,overlay.y,col="green",pch=".")
    
    title("Distribution of All Pixel values combines")
    if(ntot>maxnol){
      text(xmin+0.47*(xmax-xmin),ymin+0.9*(ymax-ymin),
           paste("Overlay of",as.character(maxnol),"of",
                 as.character(ntot),"data points",sep=" "),cex=1.3)
    }else{
      text(xmin+0.47*(xmax-xmin),ymin+0.9*(ymax-ymin),
           paste("Overlay of", as.character(ntot),"data points",sep=" "),
           cex=1.3)
    }
    text(xmin+0.47*(xmax-xmin),ymin+0.8*(ymax-ymin),
         paste("Mean =",as.character(round(mean,3)),
               "  Median =",as.character(round(median,3)),sep=" "),cex=1.3)
    text(xmin+0.47*(xmax-xmin),ymin+0.7*(ymax-ymin),
         paste("s.d. =",as.character(round(sd,3)),
               "MAD =",as.character(round(mad,3)),sep=" "),cex=1.3)
    text(xmin+0.47*(xmax-xmin),ymin+0.6*(ymax-ymin),
         paste("Gaussian(",as.character(round(median,3)),",",
               as.character(round(mad,3)),") density",sep=""),
         col="red",cex=1.3)
    if(mad>sd){
      text(xmin+0.47*(xmax-xmin),ymin+0.55*(ymax-ymin),
           paste("Warning: MAD > s.d., SigClust can be anti-conservative",
                 sep=""),cex=1.3)
    }
  }
  if(arg=="qq"|arg=="all"){
  #QQ plot
    par(mfrow=c(1,1))  
    par(mar=c(5,4,4,2)+0.1)
    qqnorm <- qqnorm(as.vector(raw.data),plot.it=FALSE)
    if(ntot>maxnol){
      which <- sample(c(1:ntot),maxnol)
    }else{
      which <- c(1:ntot)
    }
    x <- sort(qqnorm$x[which])
    y <- sort(qqnorm$y[which])
    x25 <- x[which(x>qnorm(0.25))[1]]
    x50 <- x[which(x>qnorm(0.5))[1]]
    x75 <- x[which(x>qnorm(0.75))[1]]
    y25 <- y[which(x>qnorm(0.25))[1]]
    y50 <- y[which(x>qnorm(0.5))[1]]
    y75 <- y[which(x>qnorm(0.75))[1]]
    plot(x,y,col="red",xlab="Gaussian Q",ylab="Data Q",
         main="Robust Fit Gaussian Q-Q, All Pixel values",cex.lab=1.3,...)
    abline(0,1,col="green",lwd=2)
    xmin <- min(x)-0.05*(max(x)-min(x))
    xmax <- max(x)+0.05*(max(x)-min(x))
    ymin <- min(y)-0.05*(max(y)-min(y))
    ymax <- max(y)+0.05*(max(y)-min(y))
    text(xmin+0.3*(xmax-xmin),ymin+0.9*(ymax-ymin),
         paste("Mean =",as.character(round(mean,3)),sep=" "),cex=1.3)
    text(xmin+0.3*(xmax-xmin),ymin+0.8*(ymax-ymin),
         paste("sd =",as.character(round(sd,3)),sep=" "),cex=1.3)
    text(x25,y25,"+",cex=1.3)
    text(x25+0.7,y25,"0.25 quantile",cex=1.3)
    text(x50,y50,"+",cex=1.3)
    text(x50+0.7,y50,"0.5 quantile",cex=1.3)
    text(x75,y75,"+",cex=1.3)
    text(x75+0.7,y75,"0.75 quantile",cex=1.3)
  }
  if(arg=="diag"|arg=="all"){
    #Then make Covariance Estimation Diagnostic Plot
    

    ncut <- 100
    if(d>ncut){
      par(mfrow=c(2,2))
    }else{
      par(mfrow=c(1,2))
    }
    par(mar=c(2,3.7,2,1.7))
    veigvalpos <- veigval[which(veigval > 10^(-12))]
    dpos <- length(veigvalpos)
    xmin <- 0
    xmax <- d+1
    ymin <- min(veigval) - 0.05*(max(veigval)-min(veigval))
    ymax <- max(veigval) + 0.05*(max(veigval)-min(veigval))
    plot(c(1:d),vsimeigval,type="l",lty=2,lwd=3,col="red",
         xlim=c(xmin,xmax),ylim=c(ymin,ymax),
         xlab="Component #",ylab="Eigenvalue",...)
    points(c(1:d),veigval,col="black")
    title("Eigenvalues")
    lines(c(ncut+0.5,ncut+0.5),c(ymin,ymax),col="green")

    if(icovest!=2){
      lines(c(0,d+1),c(simbackvar,simbackvar),col="magenta")
      text(xmin+0.45*(xmax-xmin),ymin+0.9*(ymax-ymin),
           paste("Background variance = ",
                 as.character(round(simbackvar,3)),sep=""),
           col="magenta")
    }
    text(xmin+0.45*(xmax-xmin),ymin+0.8*(ymax-ymin),
         "Eigenvalues for simulation",col="red")
    if(mad>sd){
      text(xmin+0.45*(xmax-xmin),ymin+0.65*(ymax-ymin),
           "Warning: MAD > s.d.",col="magenta")
    }
    
    ymin <- min(log10(veigvalpos))
    - 0.05*(max(log10(veigvalpos)) - min(log10(veigvalpos)))
    ymax <- max(log10(veigvalpos))
    + 0.05*(max(log10(veigvalpos)) - min(log10(veigvalpos)))
    plot(c(1:d),log10(vsimeigval),type="l",lty=2,lwd=3,col="red",
         xlim=c(xmin,xmax),ylim=c(ymin,ymax),
         xlab="Component #",ylab="log10(Eigenvalue)",...)
    points(c(1:dpos),log10(veigvalpos),col="black")
    title("log10 Eigenvalues")
    lines(c(ncut+0.5,ncut+0.5),c(ymin,ymax),col="green")

    if(icovest!=2){
      lines(c(0,d+1),log10(c(simbackvar,simbackvar)),col="magenta")
      text(xmin+0.45*(xmax-xmin),ymin+0.9*(ymax-ymin),
           paste("log10 Background variance = ",
                 as.character(round(log10(simbackvar),3)),sep=""),
           col="magenta")
    }
    text(xmin+0.45*(xmax-xmin),ymin+0.8*(ymax-ymin),
         "Eigenvalues for simulation",col="red")
    if(mad>sd){
      text(xmin+0.45*(xmax-xmin),ymin+0.65*(ymax-ymin),
           "SigClust may be Anti-iConservative",col="magenta")
    }
    if(length(veigval)>=ncut){
      xmin <- 0
      xmax <- ncut+1
      ymin <- min(veigval[1:ncut])
      - 0.05*(max(veigval[1:ncut]) - min(veigval[1:ncut])) 
      ymax <- max(veigval[1:ncut])
      + 0.05*(max(veigval[1:ncut]) - min(veigval[1:ncut])) 
      plot(c(1:ncut),vsimeigval[1:ncut],type="l",lty=2,lwd=3,col="red",
           xlim=c(xmin,xmax),ylim=c(ymin,ymax),
           xlab="Component #",ylab="Eigenvalue",...)
      points(c(1:ncut),veigval[1:ncut],col="black")
      title("Zoomed in version of above")

      if(icovest!=2){
        lines(c(0,d+1),c(simbackvar,simbackvar),col="magenta")
        text(xmin+0.45*(xmax-xmin),ymin+0.9*(ymax-ymin),
             paste("Background variance = ",
                   as.character(round(simbackvar,3)),sep=""),
             col="magenta")
      }
      text(xmin+0.45*(xmax-xmin),ymin+0.8*(ymax-ymin),
           "Eigenvalues for simulation",col="red")
      
      nmax <- min(dpos,ncut)
      ymin <- min(log10(veigvalpos[1:nmax])) 
      - 0.05*(max(log10(veigvalpos[1:nmax])) - min(log10(veigvalpos[1:nmax])))
      ymax <- max(log10(veigvalpos[1:nmax])) 
      + 0.05*(max(log10(veigvalpos[1:nmax])) - min(log10(veigvalpos[1:nmax])))
      plot(c(1:nmax),log10(vsimeigval[1:nmax]),type="l",lty=2,lwd=3,
           col="red",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
           xlab="Component #",ylab="log10(Eigenvalue)",...)
      points(c(1:nmax),log10(veigvalpos[1:nmax]),col="black")
      title("log10 Eigenvalues")

      if(icovest!=2){
        lines(c(0,ncut+1),log10(c(simbackvar,simbackvar)),col="magenta")
        text(xmin+0.30*(xmax-xmin),ymin+0.9*(ymax-ymin),
             paste("log10 Background variance = ",
                   as.character(round(log10(simbackvar),3)),sep=""),
             col="magenta")
      }
      text(xmin+0.30*(xmax-xmin),ymin+0.8*(ymax-ymin),
           "Eigenvalues for simulation",col="red")
    }
  }
  if(arg=="pvalue"|arg=="all"){
    
                                        # Make p Value plot
    
    par(mfrow=c(1,1))
    par(mar=c(5,4,4,2)+0.1)
    denpval <- density(simcindex)
    denrange <- quantile(denpval$y, probs=c(0,0.5, 0.75, 1))
    dy <- 0.1*(denrange[4]-denrange[1])
    
    mindex <- mean(simcindex)
    sindex <- sd(simcindex)
    
    xmin <- min(c(simcindex,xcindex))
    xmax <- max(c(simcindex,xcindex))
    dx <- 0.1*(xmax-xmin)
    xind <- seq(xmin,xmax,0.001)
    
    plot(denpval, xlim=c(xmin-dx,xmax+dx),col="red",
         xlab="Cluster Index", main="",lwd=2,...)
    title(main="SigClust Results")
    points(simcindex, runif(nsim, denrange[2], denrange[3]),
           col="blue",pch=".",cex=2)
    lines(c(xcindex, xcindex), c(denrange[1]-dy, denrange[4])+dy,
          col="green",lty=2,lwd=2)
    lines(xind,dnorm(xind,mean=mindex,sd=sindex),col="black",lty=3,lwd=2)
    legend(xmin+0.05*(xmax-xmin), denrange[4], paste("P-value=", pval),
           text.col="red",bty="n")
    legend(xmin+0.05*(xmax-xmin), denrange[3]+0.75*(denrange[4]-denrange[3]),
           paste("P-vNorm=", round(pvalnorm,3)), bty="n")
  }
}


setMethod("plot",signature(x="sigclust",y="missing"),
          function(x,y,arg="all",...){
            .plot.sigclust(x,arg,...)
          })

