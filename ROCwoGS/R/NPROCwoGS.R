NPROCwoGS <-
function (score, ncutoff, niter, CIlevel){
	
	levels(score[,3])->level
n<-length(level)
#Get the size of each Group
size<-c(1:n)
for (i in 1:n){
size[i]<-c(sum(score[,3]==level[i]))
}

#Get matrix Tscore and Rscore
matrix(0,n,max(size))->y
Tscore=y
for (i in 1:n){
	Tscore[i,c(1:size[i])]<-score[,1][score[,3]==level[i]]
	}
	
Rscore=y
for (i in 1:n){
	Rscore[i,c(1:size[i])]<-score[,2][score[,3]==level[i]]
	}
	
  G<- length(size)
  Ts <- Tscore[1,1:size[1]]
  Rs <- Rscore[1,1:size[1]]
  for(g in 2:G){
    Ts <- c(Ts,Tscore[g,1:size[g]])
    Rs <- c(Rs,Rscore[g,1:size[g]])
  }
  cutoff1<-Ts[which(Rs==1)]
  sizecutoff1<-length(cutoff1)
  cutoff2<-Ts[which(Rs==0)]
  sizecutoff2<-length(cutoff2)
  cutoff1=cutoff1%*%array(1,c(1,sizecutoff2))
  cutoff2=cutoff2%*%array(1,c(1,sizecutoff1))
  cutoff <- c(cutoff1,cutoff2)
  cutoff <- sort(cutoff)
  cutoff <- quantile(cutoff,seq(0,1,1/(ncutoff+1)))
  cutoff[c(1,ncutoff+2)]=c(min(cutoff)-1,max(cutoff)+1)
  rdirichlet <- function (n, alpha) {
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    x/as.vector(sm)
  }
  Dalpha1<- array(1/(ncutoff+1),ncutoff+1)
  Dbeta1<- array(1/(ncutoff+1),ncutoff+1)
  alpha2<-0.5
  beta2<-0.5
  theta<-0.01*array(1,G)
  pd<-array(0.05, ncutoff+1)
  DY<-array(0,c(ncutoff+1,G,2))
  U<-array(0,c(ncutoff+1,G,2))
  for(g in 1:G){
    for(i in 1:ncutoff+1){
      DY[i,g,1]<-sum(Tscore[g,1:size[g]]>=cutoff[i] & Tscore[g,1:size[g]]<cutoff[i+1] & Rscore[g,1:size[g]]==1)
      DY[i,g,2]<-sum(Tscore[g,1:size[g]]>=cutoff[i] & Tscore[g,1:size[g]]<cutoff[i+1] & Rscore[g,1:size[g]]==0)    
    }
  }
  for(g in 1:G){
    for(i in 1:ncutoff+1){
      U[i,g,1]<-rbinom(1,DY[i,g,1],theta[g]*Dbeta1[i]*(1-beta2)/(theta[g]*Dbeta1[i]*(1-beta2)+(1-theta[g])*Dalpha1[i]*alpha2))
      U[i,g,2]<-rbinom(1,DY[i,g,2],theta[g]*Dbeta1[i]*beta2/(theta[g]*Dbeta1[i]*beta2+(1-theta[g])*Dalpha1[i]*(1-alpha2)))
    }
  }
  alpha1<-1-c(0,cumsum(Dalpha1))
  beta1<- c(0,cumsum(Dbeta1))
  MC.alpha1<-  array(0,c(ncutoff+2,niter+1))
  MC.beta1<-  array(0,c(ncutoff+2,niter+1))
  MC.theta<-  array(0,c(G,niter+1))
  MC.alpha2<- array(0,niter+1)
  MC.beta2<- array(0,niter+1)
  MC.alpha1[,1]=alpha1
  MC.beta1[,1]=beta1
  MC.theta[,1]=theta
  for(iter in 1:niter){
    for(g in 1:G){
      for(i in 1:ncutoff+1){
        U[i,g,1]<-rbinom(1,DY[i,g,1],theta[g]*Dbeta1[i]*(1-beta2)/(theta[g]*Dbeta1[i]*(1-beta2)+(1-theta[g])*Dalpha1[i]*alpha2))
        U[i,g,2]<-rbinom(1,DY[i,g,2],theta[g]*Dbeta1[i]*beta2/(theta[g]*Dbeta1[i]*beta2+(1-theta[g])*Dalpha1[i]*(1-alpha2)))
      }
    }
    for(i in 1:ncutoff+1){
      pd[i]= 0.5+sum(DY[i,,]-U[i,,])
    }
    Dalpha1=  rdirichlet(1,pd)
    MC.alpha1[,iter+1]=1-c(0,cumsum(Dalpha1)[-(ncutoff+1)],1)
    for(i in 1:ncutoff+1){
      pd[i]= 0.5+sum(U[i,,])
    }
    Dbeta1=  rdirichlet(1,pd)
    MC.beta1[,iter+1]=c(0,cumsum(Dbeta1)[-(ncutoff+1)],1)
    # Gibbs Sampling alpha2
    par1<-0.5+sum(DY[,,1]-U[,,1])
    par2<-0.5+sum(DY[,,2]-U[,,2])
    alpha2<-rbeta(1, par1, par2)
    while(alpha2+beta2>1){alpha2<-rbeta(1, par1, par2)}
    MC.alpha2[iter+1]<-alpha2
    # Gibbs sampling beta2
    par1<-0.5+sum(U[,,2])
    par2<-0.5+sum(U[,,1])
    beta2<-rbeta(1, par1, par2)
    while(alpha2+beta2>1){beta2<-rbeta(1, par1, par2)}
    MC.beta2[iter+1]<-beta2
    # Gibbs Sampling theta
    for(g in 1:G){
      par1<-0.5+sum(U[,g,])
      par2<-0.5+sum(DY[,g,]-U[,g,])
      theta[g]<-rbeta(1,par1,par2)
    }
    MC.theta[,iter+1]<-theta
  }
  
  
  T.Se<-1-MC.beta1[,(niter/2+1):(niter+1)]
  T.Sp<-1-MC.alpha1[,(niter/2+1):(niter+1)]
  Prev<-MC.theta[,(niter/2+1):(niter+1)]
  R.Se<- 1-MC.beta2[(niter/2+1):(niter+1)]
  R.Sp<- 1-MC.alpha2[(niter/2+1):(niter+1)]
  flag<- 0
  if(quantile(R.Se,0.025)+quantile(R.Sp,0.025)<=1) flag=flag+1
  flag2 <- 0
  for(g in 1:(G-1)){
    for(g2 in (g+1):G){
    flag2 <- flag2 | (quantile(Prev[g,]-Prev[g2,],0.975)>=0 & quantile(Prev[g,]-Prev[g2,],0.025)<=0)
    }
  }
  if(flag2==1) flag=flag+2
  Cutoff=as.vector(cutoff[2:(ncutoff+1)]) 
  pCI<-c(0.5-CIlevel/2,    CIlevel/2+0.5)
  return(list(cutoff=Cutoff,
              T.Se=(rbind(apply(T.Se,1,mean),apply(T.Se,1,quantile, pCI))), 
              T.Sp=(rbind(apply(T.Sp,1,mean),apply(T.Sp,1,quantile, pCI))),
              R.Se=c(mean(R.Se),quantile(R.Se, pCI)),R.Sp=c(mean(R.Sp),quantile(R.Sp, pCI)),
              Prev=rbind(apply(Prev,1,mean),apply(Prev,1,quantile, pCI)),  flag=flag))
  #,quantile(Se, CIlevel/2+0.5), quantile(Se, 0.5-CIlevel/2),
  #        quantile(Sp, CIlevel/2+0.5), quantile(Se, 0.5-CIlevel/2)))
}

