#####Sub-network Selection
Subnetwork.Select=function(net,trace,node.based=NULL,infinite=TRUE,steps=5)
{
  ztall=sapply(1:ncol(trace),function(kk) return(mean(trace[,kk])))
  eids=which(ztall>0.5)
  deids=node.based
  g<-network(net)
  subnetwork=list()
  delete=NULL
  C=NULL
  
  if (length(node.based)==0){
    subnetwork$eids=eids
    subnetwork$adj=matrix(0,ncol=length(subnetwork$eids),nrow=length(subnetwork$eids))
    for (x in 1:length(subnetwork$eids))
    {
      for (y in 1:length(subnetwork$eids))
      {
        subnetwork$adj[x,y]<-net[(unlist(subnetwork$eids[x])),unlist(subnetwork$eids[y])]
      }
    }
    diag(subnetwork$adj)<-0
  }else{
    if(infinite==TRUE){
      for (i in 1:length(deids))
      {
        C=NULL
        A=get.neighborhood(g,deids[i],"combined")
        repeat{
          delete=NULL
          if (length(A)!=0){A=unique(unlist(sapply(1:length(A),function(i) get.neighborhood(g,A[i],"combined"))))}
          for (j in 1:length(A))
          {
            if (length(which(eids==A[j]))==0){delete=c(delete,j)}
          }
          if (length(delete)!=0){A=A[-delete]}
          C=unique(c(C,A))
          if (identical(sort(C),sort(Csave))) break
          Csave=C
        }
        subnetwork$eids[[i]]=unique(c(eids[i],get.neighborhood(g,deids[i],"combined"),C))
      }
    }else{
      for (i in 1:length(deids))
      {
        C=NULL
        A=get.neighborhood(g,deids[i],"combined")
        if (steps!=1){
          for(k in 1: (steps-1)){
            delete=NULL
            if (length(A)!=0){A=unique(unlist(sapply(1:length(A),function(i) get.neighborhood(g,A[i],"combined"))))}
            for (j in 1:length(A))
            {
              if (length(which(eids==A[j]))==0){delete=c(delete,j)}
            }
            if (length(delete)!=0){A=A[-delete]}
            C=unique(c(C,A))             
          }}
        
        subnetwork$eids[[i]]=unique(c(eids[i],get.neighborhood(g,deids[i],"combined"),C))
      }
    }
    
    for (i in 1:length(deids))
    {
      subnetwork$adj[[i]]=matrix(0,ncol=length(subnetwork$eids[[i]]),nrow=length(subnetwork$eids[[i]]))
      for (x in 1:length(subnetwork$eids[[i]]))
      {
        for (y in 1:length(subnetwork$eids[[i]]))
        {
          subnetwork$adj[[i]][x,y]<-net[(unlist(subnetwork$eids[[i]][x])),unlist(subnetwork$eids[[i]][y])]
        }
      }
      diag(subnetwork$adj[[i]])<-0
    }
  }
  return(subnetwork)
  
}
