pSkilMack<-function(x,b=NA,trt=NA,method=NA,n.mc=10000){
  outp<-list()
  outp$stat.name<-"Skillings-Mack SM"
  outp$n.mc<-n.mc
  
  ties<-!length(unique(as.numeric(x)))==length(x)
  
  
  if(is.matrix(x)){
    outp$n<-n<-nrow(x)
    outp$k<-k<-ncol(x)
  }
  
  ##Turn x into a matrix if it's given as a vector
  if(!is.matrix(x)){
    if ((length(x) != length(b))||(length(x) != length(trt)))
      stop("'x', 'b', and 'trt' must have the same length")
    
    outp$n<-n<-length(unique(b))
    outp$k<-k<-length(unique(trt))
    x.vec<-x
    num.obs<-length(x.vec)
    ##In case the user gives some kind of labels other than 1,2,3...
    b.ind<-as.numeric(as.factor(b))
    trt.ind<-as.numeric(as.factor(trt))
    
    ##Turn x into a matrix;
    x<-matrix(ncol=outp$k,nrow=outp$n)
    for(i in 1:num.obs){
      x[b.ind[i],trt.ind[i]]<-x.vec[i]        
    }
  }
  ###########################################
  

  outp$obs.mat<-matrix(0,ncol=outp$k,nrow=outp$n)
  outp$obs.mat[!is.na(x)]<-1
  outp$x<-x

  outp$ss<-s<-rowSums(outp$obs.mat)

  ##When the user doesn't give us any indication of which method to use, try to pick one.
  if(is.na(method)){
    if(prod(factorial(outp$ss))<=10000){
      method<-"Exact"
    }
    if(prod(factorial(outp$ss))>10000){
      method<-"Monte Carlo"
    }
  }
  #####################################################################
  
  outp$method<-method
  
  lambda.mat<-matrix(nrow=k,ncol=k)
  for(i in 1:k){
    for(j in (1:k)){
      lambda.mat[i,j]<-sum(outp$obs.mat[,i]*outp$obs.mat[,j])
    }
  }
  
  sigma.mat<-(-1)*lambda.mat[1:(k-1),1:(k-1)]
  for(i in 1:(k-1)){
    diag(sigma.mat)[i]<-sum(lambda.mat[i,-i])
  }
  
  ##Uses MASS##
  sigma0.inv<-ginv(sigma.mat)
  #############

  missing.obs<-function(rank.data){
    si<-sum(!is.na(rank.data))
    rank.data[is.na(rank.data)]<-(si+1)/2
    return(sqrt(12/(si+1))*(rank.data-(si+1)/2))
  }
  
  SM.stat<-function(obs.data){    
    ranks<-t(apply(obs.data,1,rank,na.last="keep"))
    ranks<-t(apply(ranks,1,missing.obs))
    Aj<-apply(ranks,2,sum)[1:(k-1)]
    SM.stat<-t(Aj)%*%sigma0.inv%*%t(t(Aj))
    return(as.numeric(SM.stat))
  }
  
  outp$obs.stat<-SM.stat(x);
  
  outp$obs.mat[outp$obs.mat==0]<-NA
  possible.ranks<-matrix(ncol=outp$k,nrow=outp$n)
  for(i in 1:outp$n){
    possible.ranks[i,!is.na(outp$obs.mat[i,])]<-rank(x[i,!is.na(outp$obs.mat[i,])])
  }
  
  if(outp$method=="Exact"){
    possible.perm<-multCh7SM(possible.ranks)
    exact.dist<-apply(possible.perm,3,SM.stat)
    outp$p.val<-mean(exact.dist>=outp$obs.stat)      
  }
  
  if(outp$method=="Monte Carlo"){
    mc.perm<-matrix(ncol=outp$k,nrow=outp$n)
    mc.stats<-numeric(n.mc)
    for(i in 1:n.mc){
      for(j in 1:n){
        mc.perm[j,!is.na(x[j,])]<-sample(rank(x[j,!is.na(x[j,])]))
      }
      mc.stats[i]<-SM.stat(mc.perm)
    }
    outp$p.val<-mean(mc.stats>=outp$obs.stat)      
  }  
  
  if(outp$method=="Asymptotic"){
    outp$p.val<-1-pchisq(outp$obs.stat,outp$k-1)    
  }
  class(outp)<-"NSM3Ch7p"
  outp
  
}
