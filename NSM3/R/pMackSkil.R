pMackSkil<-function(x,b=NA,trt=NA,method=NA,n.mc=10000){
  outp<-list()
  outp$stat.name<-"Mack-Skillings MS"
  outp$n.mc<-n.mc
  
  if(is.array(x)){
    outp$n<-n<-dim(x)[3]
    outp$k<-k<-dim(x)[2]
    outp$c<-dim(x)[1]
  }

  ##Turn x into an array if it's given as a vector
  if(!is.array(x)){
    if ((length(x) != length(b))||(length(x) != length(trt)))
      stop("'x', 'b', and 'trt' must have the same length")
    
    outp$n<-n<-length(unique(b))
    outp$k<-k<-length(unique(trt))
    outp$c<-length(x)/(n*k)
    
    x.vec<-x
    num.obs<-length(x.vec)
    ##In case the user gives some kind of labels other than 1,2,3...
    b.ind<-as.numeric(as.factor(b))
    trt.ind<-as.numeric(as.factor(trt))
    
    
    ##Turn x into an array;
    x<-array(dim=c(outp$c,outp$k,outp$n))
    for(i in 1:num.obs){
      tmp.ind<-1
      while(!is.na(x[tmp.ind,trt.ind[i],b.ind[i]])){tmp.ind=tmp.ind+1}
      x[tmp.ind,trt.ind[i],b.ind[i]]<-x.vec[i]        
    }
  }
  ###########################################

  num.obs<-length(x)
  
  ##When the user doesn't give us any indication of which method to use, try to pick one.
  if(is.na(method)){
    if(factorial(outp$c*outp$k)^outp$n<=10000){
      method<-"Exact"
    }
    if(factorial(outp$c*outp$k)^outp$n>10000){
      method<-"Monte Carlo"
    }
  }
  #####################################################################
  
  outp$method<-method
  
  data.mat<-matrix(x,ncol=outp$n)
  possible.ranks<-t(apply(t(data.mat),1,rank))
  
  MS.calc<-function(obs.data){
    S.vec<-NULL
    for(i in 1:outp$k){
      S.vec<-c(S.vec,sum(obs.data[,((i-1)*outp$c+1):(i*outp$c)])/outp$c)
    }
    MS.stat<-12/(outp$k*(num.obs+outp$n))*sum(S.vec^2)-3*(num.obs+outp$n)
    return(MS.stat)
  }
  
  outp$obs.stat<-MS.calc(possible.ranks)
  
  if(outp$method=="Exact"){
    possible.perm<-multCh7(possible.ranks)
    exact.dist<-apply(possible.perm,3,MS.calc)
    outp$p.val<-mean(exact.dist>=outp$obs.stat)      
  }
  if(outp$method=="Monte Carlo"){
    mc.perm<-matrix(ncol=(outp$c*outp$k),nrow=outp$n)
    mc.stats<-numeric(n.mc)
    for(i in 1:n.mc){
      for(j in 1:outp$n){
        mc.perm[j,]<-sample(possible.ranks[j,])
      }
      mc.stats[i]<-MS.calc(mc.perm)
    }
    outp$p.val<-mean(mc.stats>=outp$obs.stat)      
  }  
  
  if(outp$method=="Asymptotic"){
    outp$p.val<-1-pchisq(outp$obs.stat,outp$k-1)    
  }
  class(outp)<-"NSM3Ch7p"
  outp
}

