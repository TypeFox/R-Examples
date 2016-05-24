

MI.Cluster.Annotation=function(data,annotation,cluster){
  a.gene=annotation$GeneID
  a=annotation$Annotation
  k.gene=data$GeneID
  K=data.matrix(cluster)
  K=matrix(K,nrow=length(k.gene))
  mi=MI.score(K,k.gene,a,a.gene)
  return(mi)
}

MI.1=function(x){
  n=length(x)
  nx=table(x)+1e-20
  p=nx/n
  mi=-sum(p*log(p))
  return(mi)
}
MI.2=function(x,y){
  n=length(x)
  nxy=xtabs(~x+y)
  nx=table(x)
  ny=table(y)
  px=nx/n+1e-10
  py=ny/n+1e-10
  pxy=nxy/n+1e-10
  mi=sum(px*log(px))+sum(py*log(py))-sum(pxy*log(pxy))
  return(-mi)
}

MI.score.one=function(k,k.gene,a,a.gene){
 mi=0
 a=as.numeric(factor(a))
 for( aa in a[!duplicated(a)]){
    aa.gene=a.gene[a==aa]
    y=k.gene %in% aa.gene
    mi=mi+MI.2(k,y)
  }
 return(mi)
}

MI.score=function(K,k.gene,a,a.gene){
   if(is.vector(K)) K=matrix(K,nrow=length(K))
   bhi=rep(0,dim(K)[2])
   for( i in 1:ncol(K)){
   #  print(paste("cluster",i))
     bhi[i]=MI.score.one(K[,i],k.gene,a,a.gene)
   }
   return(bhi)
}

##---------- Normalized Mutual Information

NMI.score=function(K,k.gene,a,a.gene){
   a=as.numeric(factor(a))
   nA=length(unique(a))
   nG=length(k.gene)
   A=matrix(0,nrow=nG,ncol=nA)
   for(i in 1:nA)
     A[,i]=k.gene %in% a.gene[a==unique(a)[i]]
   
   mi=MI.score(K,k.gene,a,a.gene)
   mi.k=apply(K,2,MI.1)
   mi.a=sum(apply(A,2,MI.1))
   nmi=mi/sqrt(mi.k)/sqrt(mi.a)
   return(list(MI=mi,NMI=nmi))
}

###########################

lglk.cluster.one=function(data,model="poisson",cluster){
    n = data$Count
    s = data$Normalizer
    t = as.numeric(factor(data$Treatment))
    cluster=as.numeric(factor(cluster))

    nG=nrow(n)
    nK=max(cluster)
    nT=max(t)

   if(is.vector(s)) s=matrix(rep(s,each=nG),nrow=nG)
   if( model=="nbinom" ) nb.disp=est.nb.v.QL(n,s,t)
   if( model=="poisson") nb.disp=rep(1e-10,nG)
   if(is.vector(nb.disp)) nb.disp=matrix(rep(nb.disp,ncol(n)),nrow=nG)

   P=outer(cluster,unique(cluster),"==")+1e-10
   C=matrix(0,nrow=nK,ncol=nT)
   M=matrix(1,nrow=nG,ncol=nK)
   MC=cl.mb.est.MC(n,s,t,model,P,C,M,nb.disp)
   M=MC$M
   C=MC$C
   lglk=matrix(0,nrow=nrow(n),ncol=nK)
   for(k in 1:nK){
     r=exp(sweep(s+M[,k],2,C[k,t],"+"))
     if(model=="nbinom") 
       lglk[,k]=lglk.nb(n,r,v=nb.disp) 
     if(model=="poisson")  
       lglk[,k]=lglk.ps(n,r) 
     }
   lglk=sum(P*lglk)/nG
   return(lglk)
}
lglk.cluster=function(data,model="poisson",cluster){
   K=data.matrix(cluster)
   if(nrow(K)==1) K=t(K)
   lglk=rep(0,ncol(K))
   for(i in 1:ncol(K)){
    lglk[i]=lglk.cluster.one(data,model,K[,i])
   }
   return(lglk)
}









