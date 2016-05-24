#####DPdensity
#####DPdensity
DPdensitycluster<-function(v,rstat,DPM.mcmc,DPM.prior){
  results<-list()
  
  
  
  mcmc <- DPM.mcmc
  #prior1 <- list(alpha=0.5,m1=0.2,nu1=100,psiinv1=0.1,k0=100)
  #prior2 <- list(a0=2,b0=1,m2=rep(0,1),s2=diag(100000,1),
  #              psiinv2=solve(diag(0.5,1)),
  #             nu1=4,nu2=4,tau1=1,tau2=100)
  fit<-DPdensity(rstat,prior=DPM.prior,mcmc=mcmc,status=TRUE)
  
  for(oo in 1:v){
    nburn <-0
    nsave <-1
    nskip <-0
    ndisplay <-10
    mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)
    #prior1 <- list(alpha=0.5,m1=0.2,nu1=100,psiinv1=0.1,k0=100)
    prior2 <- list(a0=2,b0=1,m2=rep(0,1),s2=diag(100000,1),
                   psiinv2=solve(diag(0.5,1)),
                   nu1=4,nu2=4,tau1=1,tau2=100)
    fit<-DPdensity(rstat,prior=prior2,mcmc=mcmc,state=fit$state,status=FALSE)
    
    ##calculate cluster
    ss<-fit$state$ss
    num<-fit$state$ncluster
    prop<-table(fit$state$ss)/num
    mu<-fit$state$muclus[1:num]
    sigma<-fit$state$sigmaclus[1:num]
    index<-order(mu,decreasing=T)
    newm<-mu[index]
    news<-sigma[index] 
    newp<-prop[index]
    
    
    if(num==2 | num==1){results[[oo]]<-ss-1}else{
      #hh1<-1/2/sqrt(pi*news[1])+1/2/sqrt(pi*news[2])-2*dnorm(newm[1],newm[2],sd=sqrt(news[1]+news[2]))
      
      #hh2<-1/2/sqrt(pi*news[2])+1/2/sqrt(pi*news[3])-2*dnorm(newm[2],newm[3],sd=sqrt(news[2]+news[3]))
      
      clmatrix<-matrix(rep(0,num*num),ncol=num,nrow=num)
      for(i in 1:num){
        for(j in 1:num){
          clmatrix[i,j]<-dnorm(newm[i],newm[j],sd=sqrt(news[i]+news[j]))
          
        }
      }
      
      
      temp1<-0
      temp2<-0
      com1<-0
      com2<-0
      res1<-0
      res2<-0
      final1<-0
      final2<-0
      latent<-0
      latent[1]<-0
      final1<-c(1,-1,rep(0,num-2))
      final2<-c(0,-1,1,rep(0,num-3))
      
      for(iii in 1:(num-2)){
        dec1<-t(final1)%*%clmatrix%*%final1
        dec2<-t(final2)%*%clmatrix%*%final2
        
        if(dec1>dec2){latent[iii+1]<-1
        }else{
          latent[iii+1]<-0
        }
        if(latent[iii+1]==1){
          ddd<-max(which(latent==0))
          temp1<-rep(0,num)
          temp1[1:ddd]<-1
          com1<-temp1*newp
          com1<-scale(com1,center=F,scale=sum(com1))
          temp1<-rep(0,num)
          temp1[(ddd+1):(iii+1)]<--1
          res1<-temp1*newp
          res1<-scale(res1,center=F,scale=sum(abs(res1)))
          final1<-com1+res1
          temp2<-rep(0,num)
          temp2[(ddd+1):(iii+1)]<--1
          res2<-temp2*newp
          res2<-scale(res2,center=F,scale=sum(abs(res2)))
          com2<-rep(0,num)
          com2[iii+2]<-1
          final2<-com2+res2
        }else{
          temp1<-rep(0,num)
          temp1[1:iii]<-1
          com1<-temp1*newp
          com1<-scale(com1,center=F,scale=sum(com1))
          res1<-rep(0,num)
          res1[iii+1]<--1
          final1<-com1+res1
          final2<-rep(0,num)
          final2[iii+1]<--1
          final2[iii+2]<-1
        }
        
        
      }
      latent[num]<-1
      
      uuu<-max(which(latent==0))
      ss[ss %in% index[(uuu+1):num]]<-0
      ss[ss %in% index[1:uuu]]<-1
      results[[oo]]<-ss
    }
  }
  dpdensitycluster<-do.call(rbind,results)
  return(dpdensitycluster)
}
