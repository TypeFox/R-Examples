perpacf <-
function(x,T,p,missval){
       CTHRES=10^6
       ZTHRS=.Machine$double.eps*10
       nx=length(x)
         ppa=matrix(0,T,p+1)
         C=matrix(0,T,p+1)
         Splus=matrix(1,T,p+1)
         Sminus=matrix(1,T,p+1)
         nsamp=matrix(0,T,p+1)
         pmean=matrix(0,T,1)
         pmean1=matrix(0,nx,1)


         if (is.nan( missval))
            {missisnan=1
            imissx=x[(is.nan(x))]
             x[imissx]=0
            } else {
             missisnan=0
            imissx=x[x==missval]
             x[imissx]=NaN}

       for (t in 1:T)
         {index=seq(t,nx,T)
          z=x[index]
          z1=na.omit(z)
          pmean[t]=mean(z1)
          pmean1[index]=pmean[t] }

          xd=x-t(pmean1)

       for (t in 1:T)
         {baseind=seq(t,t-1,-1)
          while (min(baseind) <= 0)
         {baseind=baseind+T}

         inum=floor((nx-baseind)/T)+1
         inummin=min(inum)
         yt=xd[seq(baseind[1],nx,T)]
         yt=yt[1:inummin]

         ytminus1=xd[seq(baseind[2],nx,T)]
         ytminus1=ytminus1[1:inummin]
         nsamp[t,1]=inummin

             Rttminus1=nancorr(cbind(yt,ytminus1),1)
             Rttminus1=Rttminus1$C
             Rttminus1=as.matrix(Rttminus1)

             ppa[t,1]=Rttminus1[1,2]/sqrt(Rttminus1[2,2]*Rttminus1[1,1])
       }
             C[,1]=1

         ralpha=matrix(0,p,1)
         rbeta=matrix(0,p,1)
         for (n in 1:p)
           { for (t in 1:T)
              { baseind=seq(t+1,t-n,-1)
                while (min(baseind) <= 0)
               { baseind=baseind+T}

               inum=floor((nx-baseind)/T)+1
               inummin=min(inum)
               ytplus1=xd[seq(baseind[1],nx,T)]
               ytminus=xd[seq(baseind[n+2],nx,T)]
               ytplus1=ytplus1[1:inummin]
               ytminus=ytminus[1:inummin]
               y=matrix(0,inummin,n)

        for (s in 1:n)
            { isamp=seq(baseind[s+1],nx,T)
              isamp=isamp[1:inummin]
              y[,s]=xd[isamp]
              ralpha[s]=base::sum(na.omit(ytplus1%*%y[,s]))/inummin
              rbeta[s]=base::sum(na.omit(ytminus%*%y[,s]))/inummin
           }

        nsamp[t,n+1]=inummin

        R=nancorr(y,1)
        R=R$C
        R=as.matrix(R)

        nR=base::sum(is.nan(R[,]))
        if(nR) {R[(is.nan(R))]=0}

         C[t,n+1]=kappa(R)
         if (C[t,n+1] < CTHRES)
            { Rinv=qr.solve(R)
            }  else  {
            Rinv=gnm::MPinv(R, rank=NULL)}

        Rtp1tp1=base::sum(na.omit(ytplus1*ytplus1))/inummin
        Rtmntmn=base::sum(na.omit(ytminus*ytminus))/inummin
        Rtp1tmn=base::sum(na.omit(ytplus1*ytminus))/inummin
        signtplus1=Rtp1tp1 - ralpha[1:n]%*%Rinv%*%ralpha[1:n]
        Splus[t,n+1]=signtplus1

        signminus=Rtmntmn - rbeta[1:n]%*%Rinv%*%rbeta[1:n]
        Sminus[t,n+1]=signminus



        if(abs(signtplus1)< ZTHRS | abs(signminus) < ZTHRS)
            { ppa[t,n+1]=0
             } else {
            ppa[t,n+1]=(Rtp1tmn - t(rbeta[1:n])%*%Rinv%*%ralpha[1:n])/sqrt(signtplus1%*%signminus)
            }
    }
  }
    result = list(ppa=ppa,nsamp=nsamp)
    class(result) = "perpacf"
    result

}
