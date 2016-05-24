posterior.gamma <-
function(data,int)
  {n=length(data);thin=10;burnin=int*thin/2
   somalogdados=sum(log(data));somadados=sum(data)
   lpost=function(alpha,bbeta)
{logpost=n*(alpha*log(bbeta)-log(gamma(alpha)))+(alpha-1)*somalogdados-bbeta*somadados
logpost=logpost+(0.0001-1)*log(alpha)-0.0001*alpha+(0.0001-1)*log(bbeta)-0.0001*bbeta
       logpost}
      alphamc=array(0,c(burnin+int,1));betamc=array(0,c(burnin+int,1))
      alphamc[1]=mean(data)^2/var(data);betamc[1]=mean(data)/var(data)
      Va=(alphamc[1]/20)^2;Vb=(betamc[1]/20)^2
      for (i in 2:burnin)
         {alphaest=rgamma(1,alphamc[i-1]^2/Va,alphamc[i-1]/Va)
          betaest=rgamma(1,betamc[i-1]^2/Vb,betamc[i-1]/Vb)
          alpha=exp(lpost(alphaest,betaest)-lpost(alphamc[i-1],betamc[i-1]))
          alpha=alpha/(dgamma(alphaest,alphamc[i-1]^2/Va,alphamc[i-1]/Va))
          alpha=alpha/(dgamma(betaest,betamc[i-1]^2/Vb,betamc[i-1]/Vb))
          alpha=alpha*dgamma(alphamc[i-1],alphaest^2/Va,alphaest/Va)
          alpha=alpha*dgamma(betamc[i-1],betaest^2/Vb,betaest/Vb)
          if (is.nan(alpha)){alpha=0}
          u=runif(1)
          if (u<alpha)
            {alphamc[i]=alphaest
             betamc[i]=betaest}
          else
             {alphamc[i]=alphamc[i-1]
             betamc[i]=betamc[i-1]}
            if ((i%%100)==0)
             print(i/(burnin+thin*int))}
  alphamcb=array(0,c(int));betamcb=array(0,c(int));j=1
      for (i in (burnin+1):(burnin+thin*int))
         {alphaest=rgamma(1,alphamc[i-1]^2/Va,alphamc[i-1]/Va)
          betaest=rgamma(1,betamc[i-1]^2/Vb,betamc[i-1]/Vb)
          alpha=exp(lpost(alphaest,betaest)-lpost(alphamc[i-1],betamc[i-1]))
          alpha=alpha/(dgamma(alphaest,alphamc[i-1]^2/Va,alphamc[i-1]/Va))
          alpha=alpha/(dgamma(betaest,betamc[i-1]^2/Vb,betamc[i-1]/Vb))
          alpha=alpha*dgamma(alphamc[i-1],alphaest^2/Va,alphaest/Va)
          alpha=alpha*dgamma(betamc[i-1],betaest^2/Vb,betaest/Vb)
           if (is.nan(alpha)){alpha=0}
          if (u<alpha)
            {alphamc[i]=alphaest
             betamc[i]=betaest}
          else
             {alphamc[i]=alphamc[i-1]
             betamc[i]=betamc[i-1]}
          if ((i%%thin)==0)
              {alphamcb[j]=alphamc[i]
               betamcb[j]=betamc[i]
               j=j+1}
          if ((i%%100)==0)
             print(i/(burnin+thin*int))}
          cbind(alphamcb,betamcb)}
