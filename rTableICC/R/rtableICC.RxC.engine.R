rtableICC.RxC.engine <-
function(R,C,T,M,p,N,cluster.size,theta){
  
  if (max(cluster.size)>(length(theta))){ 
    stop(c("Maximum number of individuals in one of the clusters is ", paste(max(cluster.size)),", which is greater than maximum allowed cluster size. (1) Re-run the function,
           (2) increase maximum allowed cluster size by increasing the number of elements of theta, 
           (3) increase total number of clusters, or
           (4) decrease total number of individuals!"))    
  }
  
  rTable.raw=array(0,dim=c(R,C,N))
  dat=array(0,dim=c(max(cluster.size),R,C,M))
  rTable=array(0,dim=c(R,C))
  g.t=array(0,dim=c(R,C,(T-1)))
  g.tilde=array(0,dim=max((T-1)))
  
  pp=as.vector(t(p))  
  
  bsl=1
  say=1
  selT=0
  
  if (R==1){
    p=t(as.matrix(p))
  }else if (C==1){
    p=as.matrix(p)
  }
  sumCounts=array(0,R*C)
  for (r in 1:M){
    if (cluster.size[r]>1){
      counts=cbind(t(compositions(cluster.size[r],R*C,include.zero=TRUE)),0)
      cl=ncol(counts)
      for (j in 1:nrow(counts)){
        if (sum(counts[j,]==0)==(R*C)){
          for (i in 1:R){
            for (k in 1:C){
              if (counts[j,((i-1)*C+k)]>0){
                counts[j,cl]=theta[cluster.size[r]]*p[i,k]+(1-theta[cluster.size[r]])*p[i,k]^theta[cluster.size[r]]
              }                  
            }
          }
        } else {
          counts[j,cl]=(1-theta[cluster.size[r]])*prod(pp^counts[j,1:(R*C)])
        }
      }
      
      counts[,cl]=counts[,cl]/sum(counts[,cl])
      
      ind=rDiscrete(1,counts[,cl])$rDiscrete
      
      sel=counts[ind,1:(R*C)]
      sumCounts=sumCounts+sel
      
      if (sum(sel==0)==((R*C)-1)){
        for (i in 1:R){
          for (k in 1:C){
            if (sel[((i-1)*C+k)]>0){
              g.t[i,k,(cluster.size[r]-1)]=g.t[i,k,(cluster.size[r]-1)]+1
              rTable.raw[i,k,say:(say+sel[((i-1)*C+k)]-1)]=1
              say=say+sel[((i-1)*C+k)]
            }
          }
        }             
      } else {
        g.tilde[(cluster.size[r]-1)]=g.tilde[(cluster.size[r]-1)]+1
        
         for (i in 1:R){
           for (k in 1:C){
             if (sel[((i-1)*C+k)]>0){
               rTable.raw[i,k,say:(say+sel[((i-1)*C+k)]-1)]=1
               say=say+sel[((i-1)*C+k)]                 
             }
           }
         }
      }
    } else if (cluster.size[r]==1){
      ind=rDiscrete(1,pp)$rDiscrete
      
      sumCounts[ind]=sumCounts[ind]+1
       for (i in 1:R){
         for (k in 1:C){
           if (ind==((i-1)*C+k)){
             rTable.raw[R,C,say]=1
             say=say+1 
           }
         }
       }
    }          
  }
  rTable=sumCounts
  
  list(rTable=rTable,rTable.raw=rTable.raw,g.t=g.t,g.tilde=g.tilde)
}
