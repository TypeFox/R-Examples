tab1<-function(ratio=4,reps=100,n=1000,start=500,fbeta=FALSE,
               maxclusters=5,chains=1,clusts=c(15,15,20),
               sig=1,gamma=1,npats=80,beta=seq(-5,5,5),
               plot=FALSE){ 
  nclusts=length(clusts)
  tmp=NULL
  for (i in 1:nclusts) tmp=c(tmp,paste('beta',i,sep=""))
  pnames=c('sig2','tau2',tmp,'pcorrect')
  parm.mean=matrix(0,reps,3+nclusts)
  parm.se=matrix(0,reps,3+nclusts)
 
#   perms=permutations(nclusts)
#   nperms=dim(perms)[1]
#   tperms=t(perms)
  ngenes=sum(clusts)
  #real.perms=matrix(0,ngenes,nperms)
  mclusts=NULL
  jmatrixf<-function(jmat,maxclusters=maxclusters){
    ngenes=length(jmat)
    vs=matrix(c(jmat),ngenes,maxclusters)==
      matrix(1:maxclusters,ngenes,maxclusters,byrow=TRUE)
    return(ifelse(vs,1,0))
  }
  modef<-function(x){
    tx=table(x)
    return(as.numeric(names(tx)[which.max(tx)]))
  }
  
  for ( i in 1:nclusts){mclusts=rbind(mclusts,matrix(i,clusts[i],1))} #current labeling 1 X ngens
  realj=jmatrixf(mclusts,maxclusters=nclusts)
  diagj=diag(1/colSums(realj))
#   for (j in 1:ngenes) {for (k in 1:nperms)  real.perms[j,k]=perms[k,mclusts[j,1]]}
#   one.m=matrix(1:nclusts,nclusts,1)
#   mclusts=mclusts%*%matrix(1,1,nperms)
#   ones=matrix(1,1,nperms)
  
  for (k in 1:reps){
    ts=generate.cluster.data(ratio,npats=npats,
                             clusts=clusts,
                             sig=sig,gamma=gamma,beta=beta)
    vs=supcluster(ts,maxclusters=maxclusters,nstart=start,n=n,nchains=chains,fbeta=fbeta)
    if(plot){
      us=concordmap(vs,chains=chains)
      image(1:ngenes,1:ngenes,us$map,xlab='Genes',ylab='Genes',
                  main=paste(expression(tau^2/sigma^2),"=",ratio,"outcome included:",!fbeta),
                   col=gray(16:1 / 16))
      }
    len=n-start+1
    param=matrix(0,len*length(chains),3+nclusts)
    ms=length(vs[[1]]$parms[1,])
    for (kk in (chains)){
    for (i in 1:(n-start+1)){
      is=len*(kk-1)+i
      param[is,1]=vs[[kk]]$parms[i,1]
      param[is,2]=vs[[kk]]$parms[i,2]
      
      clustering=matrix(unlist(vs[[kk]]$parms[i,(3+maxclusters+1):ms]),1,ngenes)
      curclust=rep(0,nclusts)
      cstart=1
      for (im in 1:nclusts){
        mgenes=clusts[im]
        curclust[im]=modef(clustering[cstart:(cstart+mgenes-1)])
        cstart=cstart+mgenes
      }
      param[is,3:(2+nclusts)]=unlist(vs[[kk]]$parms[i,3+curclust]) #translating to correct beta
      correctcluster=matrix(ifelse(clustering==t(realj%*%curclust),1,0),1,ngenes)%*%realj%*%diagj
      param[is,(3+nclusts)]=mean(correctcluster)
      if (param[is,3+nclusts]<.9)  param[is,3:(2+nclusts)]=NA
      }}
    parm.mean[k,]=colMeans(param,na.rm=TRUE)
    parm.se[k,]=sqrt(diag(var(param,na.rm=TRUE)))
  }
  mean=colMeans(parm.mean)
  se=sqrt(diag(var(parm.mean)))
  mse=colMeans(parm.se)
  return(data.frame(parameter=pnames,
                    mean=mean,se=se,mse=mse))}