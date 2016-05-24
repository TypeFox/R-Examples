pDurSkiMa<-function(x,b=NA,trt=NA,method=NA,n.mc=10000){
  outp<-list()
  outp$stat.name<-"Durbin, Skillings-Mack D"
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
  
  ###Check Data structure
  if(length(unique(rowSums(!is.na(x))))!=1){
    stop("Must be same number of observations per block")
  }
  if(length(unique(colSums(!is.na(x))))!=1){
    stop("Must be same number of observations per treatment")
  }
  ###
  outp$ss<-s<-sum(!is.na(x[1,]))
  outp$pp<-p<-sum(!is.na(x[,1]))
  outp$lambda<-outp$pp*(outp$ss-1)/(outp$k-1)
  
  
  outp$obs.mat<-matrix(0,ncol=outp$k,nrow=outp$n)
  outp$obs.mat[!is.na(x)]<-1
  outp$x<-x
  
  ##When the user doesn't give us any indication of which method to use, try to pick one.
  if(is.na(method)){
    if(factorial(outp$ss)^outp$n<=10000){
      method<-"Exact"
    }
    if(factorial(outp$ss)^outp$n>10000){
      method<-"Monte Carlo"
    }
  }
  #####################################################################
  
  outp$method<-method
    
  possible.ranks<-t(apply(x,1,rank,na.last=NA))
  
  DSK.stat<-function(obs.data){
    tmp.mat<-outp$obs.mat
    for(i in 1:outp$n){
      tmp.mat[i,tmp.mat[i,]!=0]<-obs.data[i,]
    }
    Rj<-apply(tmp.mat,2,function(x) sum(x[!is.na(x)]))
    D.stat<-12/(outp$lambda*outp$k*(outp$ss+1))*sum((Rj-outp$pp*(outp$ss+1)/2)^2)
    return(D.stat)
  }
  
  outp$obs.stat<-DSK.stat(possible.ranks)
  
  if(outp$method=="Exact"){
    possible.perm<-multCh7(possible.ranks)
    exact.dist<-apply(possible.perm,3,DSK.stat)
    outp$p.val<-mean(exact.dist>=outp$obs.stat)      
  }
  if(outp$method=="Monte Carlo"){
    mc.perm<-matrix(ncol=outp$ss,nrow=outp$n)
    mc.stats<-numeric(n.mc)
    for(i in 1:n.mc){
      for(j in 1:n){
        mc.perm[j,]<-sample(possible.ranks[j,])
      }
      mc.stats[i]<-DSK.stat(mc.perm)
    }
      outp$p.val<-mean(mc.stats>=outp$obs.stat)      
  }  
  
  if(outp$method=="Asymptotic"){
    outp$p.val<-1-pchisq(outp$obs.stat,outp$k-1)    
  }
  class(outp)<-"NSM3Ch7p"
  outp
}