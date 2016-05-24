perYW <-
function(x,T,p,missval) {

  CTHRES=1e+6
  nx=length(x)
  pmean=matrix(0,T,1)
  pmean1=matrix(0,nx,1)

  ralpha<-matrix(p,1)

     if (is.nan( missval))
         {missisnan=1
          imissx=x[(is.nan(x))]
        } else {
          missisnan=0
          imissx=x[x==missval] }

          nmissx=length(imissx)

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

        pmean[i]=mean(z)            
        x[index[imiss]]=pmean[i]
        pmean1[index]=pmean[i]                                    
      }
        xd=x-pmean1


      r<-matrix(0,p,1) 
      phi<-matrix(0,T,p)
      del<-matrix(0,T,1)
      C=matrix(0,T)
      nsamp=matrix(0,T)

      for (t in 1:T)
    { baseind=seq(t,t-p,-1) 
      while (min(baseind) <= 0)         
            {baseind=baseind+T} 
       
        inum=floor((nx-baseind)/T)+1 
        inummin=min(inum)
        yt=t(xd[seq(baseind[1],nx,T)])       
        yt=(yt[1:inummin]) 
       
       y<-matrix(0,inummin,p)

      for (s in 1:p)
      { isamp=seq(baseind[s+1],nx,T)
        isamp=isamp[1:inummin]
        y[,s]=t(xd[isamp])
        ralpha[s]=sum(yt%*%y[,s]/(inummin),na.rm = TRUE)   
      }

      nsamp[t]=inummin
      nc=nancorr(y,1)
      R=nc$C
      R=as.matrix(R)
      C[t]=kappa(R) 
      if(C[t]< CTHRES)
       {Rinv=qr.solve(R) 
        } else {
        Rinv= corpcor::pseudoinverse(R)
       }
                             
      phi[t,]=ralpha%*%Rinv
      Rtt=sum(yt%*%yt)/inummin
     
      del2=Rtt-phi[t,]%*%ralpha

     if(del2>10*.Machine$double.eps)
      {del[t]=sqrt(del2) 
        } else {
       if(del2>-10*.Machine$double.eps)
       {del[t]=0 
        } else {
        stop(" del2 negative")
       }      
      }          
    }
      result = list(phi=phi, del=del)
      class(result) = "perYW"
      result
  }
