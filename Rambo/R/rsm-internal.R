.rsm.Estep <-
function(X,chi,ppi,Tau,sub,g){
  K = ncol(chi); R = max(sub); N = nrow(X); C = max(X)
  a = matrix(as.numeric(X != 0),nrow=N,byrow=F)

  digchi=digamma(chi)
  digchitot=digamma(rowSums(chi))
  digpi=digamma(ppi)
  digpitot=matrix(0,K,K)
  tz=matrix(0,N,K)
  zold=matrix(0,N,K)
  tz=Tau
  zold[,1]=zold[,1]+1
  for (k in 1:K){
    for (l in 1:K){
      digpitot[k,l]=digamma(sum(ppi[k,l,]))
    }
  }	
  
  iter=0

  while(iter<10 && sum(abs(tz-zold))>10^-3){
    iter=iter+1
    zold=tz
    for(i in 1:N){
      for(k in 1:K){
        sumterms=array(0,dim=c(N,C,K))
        for(j in 1:N){
          for(u in 1:C){
            if((a[i,j]==1 && (X[i,j]==u))||(a[j,i]==1  && (X[j,i]==u))) for (l in 1:K){
              sumterms[j,u,l]=(a[i,j]*(X[i,j]==u)*(digpi[k,l,u]-digpitot[k,l])+a[j,i]*(X[j,i]==u)*(digpi[l,k,u]-digpitot[l,k]))*tz[j,l]
            }
          }
        }
        tz[i,k]=exp(digchi[sub[i],k]-digchitot[sub[i]]+sum(sumterms))
      }
      tz[i,]=tz[i,]/sum(tz[i,])
    }
  }
  
  tz
}
.rsm.initialbis <-
function(X,sub,K){
  N = nrow(X); R = max(sub); N = nrow(X); C = max(X)
  a = matrix(as.numeric(X != 0),nrow=N,byrow=F)
  dat<-cbind(X, t(X))
  
  distance<-matrix(0, N, N)
  comneighbours<-matrix(0, N, N)
  for(i in 1:N) {
    for(j in 1:N) {
      distance[i, j]<-sum((dat[i, ] != dat[j, ]) & (dat[i, ] != 0) & (dat[j, ] != 0))
      comneighbours[i, j]<-sum((dat[i, ] != 0) & (dat[j, ] != 0))
      if(comneighbours[i, j]>0) distance[i, j]=distance[i, j]/comneighbours[i, j]
    }
  }
  diag(distance)=0
  diag(comneighbours)=0
  meanobs=sum(distance)/(sum(comneighbours>0))
  for(i in 1:N) {
    for(j in 1:N) {
      if((i != j) & (comneighbours[i, j]==0)) distance[i, j]=meanobs
    }
  }
  
  
  ztmp=array(0,dim=c(R,N,K))
  for(r in 1:R){
    
    Tau=matrix(0,N,K);
    
    rid=which(sub==r)
    ridb=which(sub!=r)
    
    assigned=seq(0,0,length=N)
    
    distpcenter=seq(0,0,length=K)
    
    meandist=sum(distance[rid,rid])/(length(rid)*(length(rid)-1))
    
    meandistinit=0
    mindistinit=0
    while(meandistinit<meandist){
      selectinit=sample(length(rid),K)
      selectinit=rid[selectinit]
      meandistinit=sum(distance[selectinit,selectinit])/(K*(K-1))
      mindistinit=min(distance[selectinit,selectinit]+diag(2*meandist,K,K))
    }
    for(k in 1:K) assigned[selectinit[k]]=k
    
    previous=seq(0,0,length=N)
    previous[rid]=-1
    cnt=0
    while(length(which(previous != assigned))>0 && cnt<10){
      cnt=cnt+1
      previous=assigned
      for (i in 1:length(rid)) if(!(rid[i] %in% selectinit)){
        assigned[rid[i]]=0
        for (k in 1:K){
          if (length(which(assigned==k))>0) distpcenter[k]=sum(distance[rid[i],which(assigned==k)])/length(which(assigned==k))
          else distpcenter[k]=0
        }
        assigned[rid[i]]=which.min(distpcenter)[1]
      }
    }
  
    previous=seq(0,0,length=N)
    cnt=0
    while(length(which(previous != assigned))>0 && cnt<10){
      cnt=cnt+1
      previous=assigned
      for(i in 1:N){
        assigned[i]=0
        for (k in 1:K){
          if (length(which(assigned==k))>0) distpcenter[k]=sum(distance[i,which(assigned==k)])/length(which(assigned==k))
          else distpcenter[k]=0
        }
        assigned[i]=which.min(distpcenter)[1]
      }
    }
    
    assignbis=matrix(0,N,2)
    for(i in 1:N) assignbis[i,]=c(i,assigned[i])
    
    Tau[assignbis]=1
    
    ztmp[r,,]=Tau;  
  }
  ztmp
}
.rsm.Mstep <-
function(X,chi,ppi,Tau,sub,g){
  K = ncol(chi); R = max(sub); N = nrow(X); C = max(X)
  a = matrix(as.numeric(X != 0),nrow=N,byrow=F)
  chi0=matrix(1,R,K); ppi0=array(1,dim=c(K,K,C));
  for(r in 1:R){
    for(k in 1:K){
      chi[r,k]=chi0[r,k]+colSums(Tau*(sub==r))[k]
    }
  }
  
  for(k in 1:K){
    for(l in 1:K){
      for(u in 1:C){
        ppi[k,l,u]=ppi0[k,l,u]+sum((X==u)*a*(t(rbind(Tau[,k]))%*%rbind(Tau[,l])))
      }
    }
  }
  
  chisum=seq(0,0,length=R)
  ppisum=matrix(0,K,K)
  
  for (r in 1:R) chisum[r]=sum(chi[r,])
  for (k in 1:K){
    for (l in 1:K) ppisum[k,l]=sum(ppi[k,l,])
  }
  
  lower=sum(lgamma(chi))-sum(lgamma(chisum))+R*lgamma(K)+sum(lgamma(ppi))-sum(lgamma(ppisum))+K*K*lgamma(C)+sum(lgamma(g))-sum(lgamma(g[,,1]+g[,,2]))-sum(.rsm.xlogx(Tau))
  
  res=list("chi"=chi, "ppi"=ppi, "lower"=lower)
  res
}
.rsm.xlogx <-
function(u){
  if(u==0 || u==1) res=0
  else res=u*log(u)
}