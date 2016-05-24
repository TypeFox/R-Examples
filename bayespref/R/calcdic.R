calcdic <-
function(lnL=NULL, popP=NULL, popM=NULL, indP=NULL, dicburn=NULL,
                  mcmc=NULL, N=NULL, pData=NULL, popN=NULL, nPops=NULL,
                  estip=NULL,measure=NULL){
  dicvec<--2 * lnL
  if (measure=="mean"){
    dbar<-mean(dicvec)
  }
  else if (measure=="median"){
    dbar<-median(dicvec)
  }
  
  if(estip==TRUE){ ## model including individual preferences
    p_cond<-0
    if (measure=="mean"){
      for (k in 1:nPops){ ## calculate conditional prob of ind preferences for each
        ## population based on mean from posterior distribution for parameters
        pp<-apply(popP[[k]][,dicburn:mcmc],1,mean)
        pm<-mean(popM[[k]][dicburn:mcmc])
        alpha<-pp * pm
        ip<-apply(indP[popN==k,,dicburn:mcmc],c(1,2),mean)
        ip<-ip/apply(ip,1,sum)
        p_cond<-p_cond + sum(log(ddirichlet(ip,alpha)))
      }
      ip<-apply(indP[,,dicburn:mcmc],c(1,2),mean)
    }
    else if (measure=="median"){
      for (k in 1:nPops){ ## calculate conditional prob of ind preferences for each
        ## population based on mean from posterior distribution for parameters
        pp<-apply(popP[[k]][,dicburn:mcmc],1,median)
        pm<-median(popM[[k]][dicburn:mcmc])
        alpha<-pp * pm
        ip<-apply(indP[popN==k,,dicburn:mcmc],c(1,2),median)
        ip<-ip/apply(ip,1,sum)
        p_cond<-p_cond + sum(log(ddirichlet(ip,alpha)))
      }
      ip<-apply(indP[,,dicburn:mcmc],c(1,2),median)
    }
    L<-0
    for(i in 1:N){ ## returns log, hence addition
      L<- L + dmultinom(pData[i,],size=sum(pData[i,]),prob=ip[i,],
                        log=TRUE)
    }
    dhat<--2 * (p_cond + L)
  }

  else { ## model excluding individual preferences
    L<-0
    if (measure=="mean"){
      for(i in 1:N){
        curpop<-popN[i]
        pp<-apply(popP[[curpop]][,dicburn:mcmc],1,mean)
        pm<-mean(popM[[curpop]][dicburn:mcmc])
        alpha<-pp * pm
        L<- L + log(dpolya(pData[i,],alpha=alpha))
      }
    }
    else if (measure=="median"){
      for(i in 1:N){
        curpop<-popN[i]
        pp<-apply(popP[[curpop]][,dicburn:mcmc],1,median)
        pm<-median(popM[[curpop]][dicburn:mcmc])
        alpha<-pp * pm
        L<- L + log(dpolya(pData[i,],alpha=alpha))
      }
    }
    dhat<--2 * L
  }

  
  pD<-dbar - dhat
  dic<-pD + dbar
  return(dic)
}

