merge.chains=function(x,chains=NULL){
  if (is.null(chains)) chains=(1:length(x))
  outpt=NULL
  for (i in chains) outpt=rbind(outpt,x[[i]]$parms)
  return(list(inp=x[[1]]$inp,parms=outpt))}

beta.by.gene=function(supcluster.list){
  x=supcluster.list
  m=length(x)
  maxclusters=x[[1]]$inp[1]
  ngenes=x[[1]]$inp[2]
  vs1=NULL
  for (i in 1:m){
    nn=dim(x[[1]]$parms)[1]
    vs=matrix(0,nn,ngenes)
 for(j in (1:nn)) vs[j,]=unlist(x[[i]]$parms[j,3+unlist(x[[i]]$parms[j,(3+maxclusters+1):(3+maxclusters+ngenes)]) ])
 vs1=rbind(vs1,cbind(rep(i,nn),vs))
  }
 return(vs1)
 }

concordmap<-function(supcluster.list,chains=1,sort.genes=FALSE,criteria=1){
  jmatrixf<-function(jmat,maxclusters){
    ngenes=length(jmat)
    vs=matrix(c(jmat),ngenes,maxclusters)==
      matrix(1:maxclusters,ngenes,maxclusters,byrow=TRUE)
    return(ifelse(vs,1,0))
  }
  otpt=merge.chains(supcluster.list,chains=chains)
  dimotpt=dim(otpt$parms)
  maxclusters=otpt$inp[1]
  ngenes=otpt$inp[2]
  ms=dim(otpt$parms)[1]
  sameclust=matrix(0,ngenes,ngenes)
  for (j in (1:ms)){
    jmatrix=jmatrixf(otpt$parms[j,(3+maxclusters+1):(3+maxclusters+ngenes)],maxclusters)
    sameclust=sameclust+tcrossprod(jmatrix,jmatrix)
  }
  ordering=1:ngenes
  sameclust=sameclust/ms
  clusters=NULL
  if (sort.genes){
    ordering=NULL
    sameclust1=sameclust-diag(1,ngenes,ngenes)
    #genes that are in a cluster
    inclusts=(1:ngenes)[(colSums(sameclust1>=criteria)>0)]
    ninclusts=setdiff(1:ngenes,inclusts)
    nclusts=NULL
    j=1
    while((length(order)<ngenes)&&length(inclusts>0)){
      frst=min(inclusts)
      nxt=(1:ngenes)[sameclust[frst,]==1]
      ordering=c(ordering,nxt)
      clusters=c(clusters,rep(otpt$parms[ms,3+maxclusters+frst],length(nxt)))
      nclusts=(c(nclusts,rep(j,length(nxt))))
      j+1
      inclusts=setdiff(inclusts,ordering)    
    }
    lclust=length(clusters)
    if (lclust<ngenes){
      ranking=sameclust[ninclusts,ordering]%*%matrix(nclusts,length(ordering),1)
      ordering=c(ordering,ninclusts[order(ranking)])
      clusters=c(clusters,rep(NA,ngenes-length(clusters)))
    }
  } else ordering=1:ngenes
  return(list(map=sameclust[ordering,ordering],order=ordering,clusters=clusters))
}

compare.chains<-function(supcluster.list,chains1,chains2){
  map1=concordmap(supcluster.list,chains=chains1)$map
  map2=concordmap(supcluster.list,chains=chains2)$map
  ms=dim(map1)[1]
  pairs=NULL
  for (i in 1:(ms-1)){for (j in (i+1):(ms-1)){
    pairs=rbind(pairs,c(i,j,map1[i,j],map2[i,j]))
  }
  }
  return(pairs)
}







