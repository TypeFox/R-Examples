posterior.gumbel <-
function(data,block,int)
  {thin=10;burnin=int*thin/2  
   ajuste=gev(data,block)
   data=ajuste$data;n=length(data)
   lpost=function(mu,sigma)
      {logpost=-n*log(sigma)-sum(((data-mu)/sigma))-sum(exp(-(data-mu)/sigma))
       logpost=logpost+(0.001-1)*log(sigma)-0.001*sigma-mu^2/2000
       logpost}
      mumc=array(0,c(burnin+int,1));sigmamc=array(0,c(burnin+int,1))
      mumc[1]=ajuste$par.ests[3];sigmamc[1]=ajuste$par.ests[2]
      Vu=(sigmamc[1]/10)
      Vsigma=(sigmamc[1]/25)^2
      for (i in 2:burnin)
         {muest=rnorm(1,mumc[i-1],Vu)
          sigmaest=rgamma(1,sigmamc[i-1]^2/Vsigma,sigmamc[i-1]/Vsigma)
      alpha=exp(lpost(muest,sigmaest)-lpost(mumc[i-1],sigmamc[i-1]))
      alpha=alpha*dgamma(sigmamc[i-1],sigmaest^2/Vsigma,sigmaest/Vsigma)
      alpha=alpha/(dgamma(sigmaest,sigmamc[i-1]^2/Vsigma,sigmamc[i-1]/Vsigma))
       if (is.nan(alpha)){alpha=0}  
          u=runif(1)
          if (u<alpha)
            {mumc[i]=muest
             sigmamc[i]=sigmaest}
          else
             {mumc[i]=mumc[i-1]
             sigmamc[i]=sigmamc[i-1]}
            if ((i%%100)==0)
             print(i/(burnin+thin*int))}
  mumcb=array(0,c(int));sigmamcb=array(0,c(int));j=1
      for (i in (burnin+1):(burnin+thin*int))
         {muest=rnorm(1,mumc[i-1],Vu)
          sigmaest=rgamma(1,sigmamc[i-1]^2/Vsigma,sigmamc[i-1]/Vsigma)
      alpha=exp(lpost(muest,sigmaest)-lpost(mumc[i-1],sigmamc[i-1]))
      alpha=alpha*dgamma(sigmamc[i-1],sigmaest^2/Vsigma,sigmaest/Vsigma)
      alpha=alpha/(dgamma(sigmaest,sigmamc[i-1]^2/Vsigma,sigmamc[i-1]/Vsigma))
       if (is.nan(alpha)){alpha=0}
          u=runif(1)
          if (u<alpha)
            {mumc[i]=muest
             sigmamc[i]=sigmaest}
          else
             {mumc[i]=mumc[i-1]
             sigmamc[i]=sigmamc[i-1]}
          if ((i%%thin)==0)
              {mumcb[j]=mumc[i]
               sigmamcb[j]=sigmamc[i]
               j=j+1}
          if ((i%%100)==0)
             print(i/(burnin+thin*int))}
          cbind(mumcb,sigmamcb)}
