ppfplot <-
function(ppf,nsamp,alpha,datastr){

     T=nrow(ppf)
     nc=ncol(ppf)
     thr=qnorm(1-alpha/2,0,1/sqrt(nsamp))
     bfalpha=alpha/(T*nc)
     bfthr=qnorm(1-bfalpha/2,0,1/sqrt(nsamp))
     ylab.name=expression( pi(t,n+1))


      dev.set(which=1)
      ylab.name=expression( pi(t,n+1))
      matplot(seq(0,nc-1),t(ppf), xlab="n = samples between", ylab=ylab.name,type="l",lwd=1, ylim=c(-1,1))

      lines(seq(0,nc-1),-thr%*%matrix(1,1,nc),type="l", col="red",lwd=1)
      lines(seq(0,nc-1),thr%*%matrix(1,1,nc),type="l",  col="red",lwd=1)
      lines(seq(0,nc-1),-bfthr%*%matrix(1,1,nc),type="l",  col="black",lwd=1)
      lines(seq(0,nc-1),bfthr%*%matrix(1,1,nc),type="l", col="black",lwd=1)

      ccs<-matrix(0,1,T)
      for (i in 1:T)  { ccs[i]=i}
      ccs=as.character(ccs)

      par(xpd=NA,oma=c(3,0,0,0))
      legend(nc-2,1.5,ccs, col=ccs, lty=1,title="t = seasons")
      title(paste('perpacf for ', datastr))


     dev.set(which=1)
     par(mfrow = c(4,ceiling(T/4)))
     for (i in 1:T)
    { plot(seq(0,nc-1),ppf[i,], xlab="n = samples between ", t = "h", ylab=ylab.name,lwd=1, ylim=c(-1,1))
      abline(h = 0, col = "black")
      lines(seq(0,nc-1),-thr%*%matrix(1,1,nc),type="l", col="red",lwd=1)
      lines(seq(0,nc-1),thr%*%matrix(1,1,nc),type="l",  col="red",lwd=1)
      lines(seq(0,nc-1),-bfthr%*%matrix(1,1,nc),type="l",  col="blue",lwd=1)
      lines(seq(0,nc-1),bfthr%*%matrix(1,1,nc),type="l", col="blue",lwd=1)
      title(paste( 'season =',i ))  }


     colmeans<-matrix(0,1,nc)
     for (i in 1:nc)  { colmeans[i]=mean(ppf[,i]) }

     colmeans<-matrix(0,1,nc)
     for (i in 1:nc)  { colmeans[i]=mean(ppf[,i]) }

      cat(paste('p-values for pi(t,n+1)=0 for t=0,...,T-1','\n'))
      cat(paste(' for nmax=',nc,'\n'))
      nrows=nc
      pv1save<-0

      pv1<-matrix(1,nc,0)
      pv2<-matrix(1,nc,0)

      for (j in 1:nc)
        {  mean_T=colmeans[j]
           mzscore1=mean_T%*%sqrt(nsamp)%*%sqrt(T)
           pv1[j]=2*(1-pnorm(abs(mzscore1),0,1))
           pv1save=cbind(pv1save,pv1[j])
           mzscore2=mean(colmeans[j:nc])%*%sqrt(nsamp)%*%sqrt(T)%*%sqrt(nrows)
           pv2[j]=2*(1-pnorm(abs(mzscore2),0,1))
           if (j==1) {pv2save=pv2[j]}
           nrows=nrows-1
        }


            detail <- matrix(c(pv1,pv2),ncol=2)
            colnames(detail) <- c( "pi(:,n+1)=0", "pi(:,n+1:nmax)=0")
            row.names(detail)<-paste("n=",seq(0,nc-1), sep="")
             print(detail)


       pv1save=pv1save[2:length(pv1save)]
       cat(paste('test for all pi(t,n+1)=0 (1st row, 2nd col) pv=',pv2save,'\n'))
       cat(paste('min of first col, bf corrected for',nc,' rows:',nc*min(pv1save),'\n'))
 }
