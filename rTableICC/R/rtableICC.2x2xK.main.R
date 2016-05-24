rtableICC.2x2xK.main <-
function(p,theta,M,sampling="Multinomial",N=0,lambda=2,zero.clusters=FALSE,print.regular,print.raw){
  K=nrow(p)
  num.cell=2*2*K
  num.centers=length(M)                     
  if (num.centers==1){
    M=rep(M,K)
    num.centers=K 
  } else if (num.centers!=K){
    stop("Length of M should be either 1 or K!")
  }
  num.cluster=M
  if (sampling=="Product"){
    if ((length(N)!=K) | (is.finite(N)==FALSE)){ 
      stop("Total number of observations in each center should be entered under product-multinomial samlping plan. Please enter a Kx1 dimensional vector as N!")      
    }else{
      N=abs(round(N))  
      rTable.raw=array(0,dim=c(sum(N),2*K,2))
      total.N=sum(N)
      
      p=sweep(p, 2, (N/total.N), "/")
      
      scaled.p=p/apply(p,1,sum)     
      
      cluster.size=array(0,dim=c(num.centers,max(num.cluster)))
      if (zero.clusters==TRUE){
        for (i in 1:num.centers){ 
          cluster.size[i,1:num.cluster[i]]=rmultinom(1, N[i], rep(1/num.cluster[i],num.cluster[i]))        
        }
      }else if (zero.clusters==FALSE){
        for (i in 1:num.centers){ 
          for (j in 1:N[i]){
            if ((N[i]-num.cluster[i])<0){
              stop("Because number of individuals is less than the total number of clusters, it is impossible to allocate an individual to each cluster! Set zero.clusters = TRUE and re-run the function.")
            }
            cluster.size[i,1:num.cluster[i]]=rmultinom(1, (N[i]-num.cluster[i]), rep(1/num.cluster[i],num.cluster[i]))+1
          }
        }
      }
    }    
  }else if (sampling=="Multinomial"){      
    if ((length(N)!=1) | (is.finite(N)==FALSE)){ 
      stop("Total number of observations must be entered as a scalar greater than zero under multinomial samlping plan!")      
    }else{
      N=abs(round(N))
      rTable.raw=array(0,dim=c(N,2*K,2))
      cluster.size=0
      if (zero.clusters==TRUE){
        cluster.size=rmultinom(1, N, rep(1/sum(num.cluster),sum(num.cluster)))
      }else if (zero.clusters==FALSE){
        if ((N-sum(num.cluster))<0){
          stop("Because number of individuals is less than the total number of clusters, it is impossible to allocate an individual to each cluster! Set zero.clusters = TRUE and re-run the function.")
        }
        cluster.size=rmultinom(1, (N-sum(num.cluster)), rep(1/sum(num.cluster),sum(num.cluster)))+1
      }

      scaled.p=p/apply(p,1,sum)  
      N=rep(N,num.centers)
    }
  }else if (sampling=="Poisson"){
    scaled.p=p/apply(p,1,sum) 
    cluster.size=0
    if (length(lambda)>1){
      if (zero.clusters==TRUE){
        cluster.size=rpois(num.cluster[1],lambda[1])
        for (i in 2:num.centers){
          cluster.size=c(cluster.size,rpois(num.cluster[i],lambda[i]))
        }
      }else if (zero.clusters==FALSE){
        cluster.size=rktp(num.cluster[1],0,lambda[1],xpred = 1)
        for (i in 2:num.centers){
          cluster.size=c(cluster.size,rktp(num.cluster[i],0,lambda[i],xpred = 1))
        }
      }
    }else{
      if (zero.clusters==TRUE){
        cluster.size=rpois(sum(num.cluster),lambda)
      }else if (zero.clusters==FALSE){
        cluster.size=rktp(sum(num.cluster),0,lambda,xpred = 1)
      }
    }
    rTable.raw=array(0,dim=c(sum(cluster.size),2*K,2))
    N=sum(cluster.size)
  }
  if (max(cluster.size)>(length(theta))){
    stop(c("Maximum number of individuals in one of the clusters is ", paste(max(cluster.size)),", which is greater than maximum allowed cluster size. (1) Re-run the function,
           (2) increase maximum allowed cluster size by increasing the number of elements of theta, 
           (3) increase total number of clusters, or
           (4) decrease total number of individuals!"))
  }
  dat=array(0,dim=c(max(num.cluster),6,num.centers))
  rTable=array(0,dim=c(K,4))
  g.t=array(0,dim=c((2*K),2,max(N)-1))
  g.tilde=array(0,dim=max(N))
  
  bsl=1
  say=1
  selT=0
  for (i in 1:num.centers){
    if (sampling=="Product"){
      dat[1:num.cluster[i],1,i]=cluster.size[i,1:num.cluster[i]]
    }else if ((sampling=="Multinomial") | (sampling=="Poisson")){
      dat[1:num.cluster[i],1,i]=cluster.size[bsl:(bsl+num.cluster[i]-1)]      
      bsl=bsl+num.cluster[i]
    }   
    
    sumCounts=array(0,4)
    for (r in 1:num.cluster[i]){
      if (dat[r,1,i]>1){ 
        counts=0
        counts=t(compositions(dat[r,1,i],4,include.zero=TRUE)) 
        counts=cbind(counts,array(0,dim=c(nrow(counts),1)))
        
        for (j in 1:nrow(counts)){
          if (sum(counts[j,]==0)==4){
            if (counts[j,1]>0){
              counts[j,5]=theta[dat[r,1,i]]*scaled.p[i,1]+(1-theta[dat[r,1,i]])*scaled.p[i,1]^theta[dat[r,1,i]]
            }else if (counts[j,2]>0){
              counts[j,5]=theta[dat[r,1,i]]*scaled.p[i,2]+(1-theta[dat[r,1,i]])*scaled.p[i,2]^theta[dat[r,1,i]]
            }else if (counts[j,3]>0){
              counts[j,5]=theta[dat[r,1,i]]*scaled.p[i,3]+(1-theta[dat[r,1,i]])*scaled.p[i,3]^theta[dat[r,1,i]]
            }else if (counts[j,4]>0){
              counts[j,5]=theta[dat[r,1,i]]*scaled.p[i,4]+(1-theta[dat[r,1,i]])*scaled.p[i,4]^theta[dat[r,1,i]]
            }
          } else {
            counts[j,5]=(1-theta[dat[r,1,i]])*prod(scaled.p[i,]^counts[j,1:4])
          }
        }
        counts=cbind(counts,r,i)
      
        counts[,5]=counts[,5]/sum(counts[,5])
 
        ind=rDiscrete(1,counts[,5])
        ind=ind$rDiscrete
        sel=counts[ind,1:4]
        sumCounts=sumCounts+sel
        if (sum(sel==0)==3){
          if (sel[1]>0){
            g.t[2*i-2+1,1,dat[r,1,i]-1]=g.t[2*i-2+1,1,dat[r,1,i]-1]+1  
            rTable.raw[say:(say+sel[1]-1),2*i-2+1,1]=1
            say=say+sel[1]
          }else if (sel[2]>0){
            g.t[2*i-2+1,2,dat[r,1,i]-1]=g.t[2*i-2+1,2,dat[r,1,i]-1]+1
            rTable.raw[say:(say+sel[2]-1),2*i-2+1,2]=1
            say=say+sel[2]
          }else if (sel[3]>0){
            g.t[2*i-2+2,1,dat[r,1,i]-1]=g.t[2*i-2+2,1,dat[r,1,i]-1]+1
            rTable.raw[say:(say+sel[3]-1),2*i-2+2,1]=1
            say=say+sel[3]
          }else if (sel[4]>0){
            g.t[2*i-2+2,2,dat[r,1,i]-1]=g.t[2*i-2+2,2,dat[r,1,i]-1]+1
            rTable.raw[say:(say+sel[4]-1),2*i-2+2,2]=1
            say=say+sel[4]
          }              
        } else {
          g.tilde[dat[r,1,i]]=g.tilde[dat[r,1,i]]+1
          if (sel[1]>0){
            rTable.raw[say:(say+sel[1]-1),2*i-2+1,1]=1
            say=say+sel[1]                 
          }
          if (sel[2]>0){
            rTable.raw[say:(say+sel[2]-1),2*i-2+1,2]=1
            say=say+sel[2]                 
          }
          if (sel[3]>0){
            rTable.raw[say:(say+sel[3]-1),2*i-2+2,1]=1
            say=say+sel[3]                 
          }
          if (sel[4]>0){
            rTable.raw[say:(say+sel[4]-1),2*i-2+2,2]=1
            say=say+sel[4]                 
          }
        }
      } else if (dat[r,1,i]==1){
        ind=rDiscrete(1,scaled.p[i,])
        ind=ind$rDiscrete
        sumCounts[ind]=sumCounts[ind]+1
        if (ind==1){
          rTable.raw[say,2*i-2+1,1]=1
          say=say+1  
        } else if (ind==2){
          rTable.raw[say,2*i-2+1,2]=1
          say=say+1  
        } else if (ind==3){
          rTable.raw[say,2*i-2+2,1]=1
          say=say+1  
        } else if (ind==4){
          rTable.raw[say,2*i-2+2,2]=1
          say=say+1  
        }
      }          
    }    
    rTable[i,]=sumCounts
  }
  T=max(cluster.size)
  g.t=g.t[,,1:(T-1)]
  g.tilde=g.tilde[1:(T-1)]
  if (print.regular==TRUE){
    rTable.regular=array(0,dim=c(2,2,K))
    for (i in 1:K){
      say=0
      for (j in 1:2){
        for (r in 1:2){
          say=say+1
          rTable.regular[j,r,i]=rTable[i,say]      
        }
      }
    }
    list(g.t=g.t,g.tilde=g.tilde,rTable=rTable,rTable.raw=rTable.raw,rTable.regular=rTable.regular,N=N,cluster.size=cluster.size,sampling=sampling,
         M=M,K=K,T=T,ICC=TRUE,structure="2x2xK",print.raw=print.raw,print.regular=print.regular)
  }else {  
     list(g.t=g.t,g.tilde=g.tilde,rTable=rTable,rTable.raw=rTable.raw,N=N,cluster.size=cluster.size,sampling=sampling,
       M=M,K=K,T=T,ICC=TRUE,structure="2x2xK",print.raw=print.raw,print.regular=print.regular)
  }
}
