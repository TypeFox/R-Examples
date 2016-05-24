iFCB<- function(data,nclus=3,ndim=2,nstart=100,smartStart=F,seed=1234){
  data = data.matrix(data) 
  outDisj=disjMake(data)
  
  #in case of binary input
  minobs = min(data)
  maxobs = max(data)
  
  q=ncol(data)
  n=nrow(data)
  dZ=outDisj$dZ
  
  Q=ncol(dZ)  
  best_f=1000000
  
  for(b in 1:nstart){
    
    set.seed(seed+b)
    
    gvec =matrix(ceiling(runif(n)*nclus),n,1)
    #print(ceiling(runif(n)*nclus)[1:10])
    if(smartStart==T){
      outkstart=kmeans(dZ,nclus,nstart=100)
      gvec=outkstart$cluster  
    }
    C=matrix(0,n,nclus)
    for (i in 1:n){
      C[i,gvec[i]]=1
    }
    
    #print("random cluster allocation matrix")
    w= -Inf
    ceps=0.00000001
    itmax=100
    it=0 ### inizializzazione contatore iterazioni
    imp=100000 
    f0 =10000000 ## inizializzazione criterio di arresto
    u=matrix(1,n,1);## vettore di 'uno' 
    #story.obscoord=list()
    
    while((it<=itmax)&&(imp>ceps)){
      
      it=it+1 
      #print("it")
      #print(it)
      
      Fmat=t(C) %*% dZ
      P=Fmat/sum(Fmat)
      r= apply(P,1,sum)
      c= apply(P,2,sum)
      
      r=t(t(r))
      c=t(t(c))
      onec=matrix(1,nrow=ncol(P))
      nsSpc=sqrt(q)*   t(t(t(t(P)*as.vector(1/c)) - r %*% t(onec)) * as.vector(sqrt(c)))
   #
      
      nssvdres=svd(nsSpc)
      
     # print("decomposition done")
      
      nU  = nssvdres$u[,1:ndim]
      nV  = nssvdres$v[,1:ndim]
      nsv = nssvdres$d[1:ndim]
      
      G =  t(t(nU)*nsv)
      
      B =   (1/(sqrt(as.vector(c)*sqrt(q)))) * t(t(nV) * nsv)  
      
      Csize = apply(C,2,sum)
      Cw = as.vector(C %*% Csize)
      Y = Cw * (dZ / (n*sqrt(q))) %*% t(t(B) * as.vector(1/nsv))
      
      Kout = kmeans(Y,centers=G)
      G=Kout$centers
      ngvec = Kout$cluster
      
      C=matrix(0,n,nclus)
      for (i in 1:n){
        C[i,ngvec[i]]=1
      }
  
      centerC=matrix(apply(C,2,sum),nrow=n,ncol=nclus,byrow=T)/n
      centerZ=matrix(c*(n*q),nrow=n,ncol=ncol(dZ),byrow=T)/n
  
      Cstar=C-centerC
      Zstar=dZ-centerZ
      
      fA=Cstar- t(t(Zstar) * as.vector(sqrt(c))) %*% B %*% t(nU) 
      flossA=sum(diag(t(fA) %*% fA))  
      flossB=sum(diag(t(Y - Cstar %*% G) %*% (Y - Cstar %*% G)))    
      f=flossA+flossB 
      imp=f0-f
      f0=f
    }
    
    if(f<=best_f){
      
      best_f=f
      best_ngvec=ngvec	
      best_Y=Y
      best_B=B
      best_G=G
      best_it=it
    }
    
    
    
  }## end FOR
  
  #####################
  ####################
  f=best_f
  ngvec=best_ngvec	
  Y=best_Y	
  B=best_B
  G=best_G
  
  it=best_it
  
  ####################
  #####################
  if ((minobs==0) & (maxobs==1)) {
    B=B[seq(from=1,to=(2*q),by=2),]
  }
  otpt=list() 
  otpt$obscoord=Y
  otpt$attcoord=B
  otpt$centroid=G
  otpt$cluID=ngvec
  otpt$criterion=f
  otpt
}