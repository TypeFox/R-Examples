permest <-function(x,T,alpha,missval,datastr,...){

 permest_full<-function(x,T,alpha,missval,datastr,typeci,typepmean,pchci,pchpmean,colci,colpmean,pp)
 {
       nx=length(x)
       nper = floor (nx/T)
       nxact=nper*T
       nrem=nx-nxact
       if ( nrem>0) {nper=nper+1}


       if (is.nan( missval)) {
          missisnan=1
          imissx=x[(is.nan(x))]
         } else {
          missisnan=0
          imissx=x[x==missval] }

      nmissx=length(imissx)


      pmean1<-matrix(0,nx,1)
      pmean<-matrix(0,T,1)
      pstd<-matrix(0,T,1)
      ny<-matrix(0,T,1)
      pmci<-matrix(0,T,2)

      bigz=c()
      groupz=c()
      X=NaN*matrix(1,nper,T)

    if (pp)
     {
      cat(paste('found ',nper,' periods of length ',T,' with remainder of ',nrem,'\n'))
     }

     vimiss<-c()
     for (i in 1:T)
     {  index=seq(i,nx,T)
               z=x[index]
               if (missisnan)
              { igood=which(!is.nan(z))
                imiss=which(is.nan(z))
                }  else  {
                igood=which((z!=missval))
                imiss=which((z==missval))
              }

         z=z[igood]
         bigz<-c(bigz,z)

        ny[i]=length(z)
        if (ny[i]==nper)
           {X[,i]=z
            } else  {
           X[(1:ny[i]),i]=z
           X[(ny[i]+1):nper,i]=NaN }

        if ( ny[i]>0)
         { groupz<-c(groupz,i*matrix(1,ny[i],1)) }

      pmean[i]=mean(z)
      pstd[i]=sd(z)
      x[index[imiss]]=pmean[i]
      xr=t(x)
      pmean1[index]=pmean[i]

       t0 = qt(c(alpha/2, 1-alpha/2),ny[i]-1)
       pmci[i,1] = pmean[i] + t0[1]*pstd[i]/sqrt(ny[i]-1)
       pmci[i,2] = pmean[i] + t0[2]*pstd[i]/sqrt(ny[i]-1)

      vimiss[i]=length(imiss)
}

  if (pp)
     { detail <- matrix(c(ny,vimiss,pmean,pmci[,1], pmci[,2]),ncol=5)
       colnames(detail) <- c(" ngood", "nmiss"," pmean ", "lower", "upper")
       row.names(detail)<-paste("i=",seq(1,T), sep="")
       print(detail) }

    xd=t(x)-pmean1

    pmean.aov= aov(bigz~groupz)
    s<-summary(pmean.aov)
    st=t(s[[1]])
    pmpv=st[5]

   if (pp)
       { matplot(pmci, xlab="seasons", ylab="mean", type=typeci,lwd=1, lty=1, col=colci,pch=pchci)
         points(pmean,type=typepmean,lwd=1, lty=1, col=colpmean,pch=pchpmean)
         title(main=(paste("Periodic mean: ","No. periods =", nper," alpha =", alpha)),sub=(paste("anova p-value for m(t) = m:",pmpv)))
         legend("bottomright", c(expression(mean),expression(confidence_intervals)), fill=c(colpmean,colci),ncol=2,title="legend")
       }
    result = list(xr=xr, xd= xd, pmean = pmean, pmci = pmci, pmpv = pmpv)
    class(result) = "permest"
    result
}

L<-modifyList(list(typeci="o",typepmean="b",pchci=10,pchpmean=15,colci="red",colpmean="blue",pp=1), list(x = x, T=T, alpha=alpha, missval=missval, datastr=datastr,...))

 do.call(permest_full,L)

 }
