pNWWM<-function(x,b=NA,trt=NA,method=NA, n.mc=10000){
  outp<-list()
  outp$stat.name<-"Nemenyi, Wilcoxon-Wilcox, Miller R*"
  outp$n.mc<-n.mc
  
  ties<-!length(unique(as.numeric(x)))==length(x)
  
  if(is.matrix(x)){
    outp$n<-n<-nrow(x)
    outp$k<-k<-ncol(x)
  }
  if(!is.matrix(x)){
    if ((length(x) != length(b))||(length(x) != length(trt)))
      stop("'x', 'b', and 'trt' must have the same length")
    
    outp$n<-n<-length(unique(b))
    outp$k<-k<-length(unique(trt))
    x.vec<-x
    ##In case the user gives some kind of labels other than 1,2,3...
    b.ind<-as.numeric(as.factor(b))
    trt.ind<-as.numeric(as.factor(trt))
    ##Turn x into a matrix;
    x<-matrix(ncol=outp$k,nrow=outp$n)
    for(i in 1:outp$n){
      for(j in 1:outp$k){
        x[i,j]<-x.vec[(b==i)&(trt==j)]        
      }
    }
  }
  
  ##When the user doesn't give us any indication of which method to use, try to pick one.
  if(is.na(method)){
    if(outp$k*factorial(outp$k)^outp$n<=10000){
      method<-"Exact"
    }
    if(outp$k*factorial(outp$k)^outp$n>10000){
      method<-"Monte Carlo"
    }
  }
  #####################################################################
  
  outp$method<-method  
  
  R.star.calc<-function(two.dim.mat,u){
    row.ranks<-t(apply(two.dim.mat,1,rank))
    return(colSums(row.ranks)[u]-colSums(row.ranks)[1])
  }
  
  R.star.all<-function(two.dim.mat){
    row.ranks<-t(apply(two.dim.mat,1,rank))
    return(max(colSums(row.ranks)[-1])-colSums(row.ranks)[1])      
  }
  
  outp$num.comp<-num.comp<-outp$k-1 
  
  count<-1
  outp$labels<-character(outp$num.comp)
  outp$obs.stat<-outp$p.val<-numeric(outp$num.comp)
  
  for(j in 2:outp$k){
    outp$labels[count]<-paste("1-",j)
    outp$obs.stat[count]<-R.star.calc(x,j)
    count<-count+1
  }  
  
  possible.ranks<-t(apply(x,1,rank))
  
  if(outp$method=="Exact"){
    possible.perm<-multCh7(possible.ranks)
    exact.dist<-apply(possible.perm,3,R.star.all)
    for(i in 1:outp$num.comp){
      outp$p.val[i]<-mean(exact.dist>=outp$obs.stat[i])      
    }
  }
  
  
  if(outp$method=="Monte Carlo"){
    mc.perm<-matrix(ncol=outp$k,nrow=outp$n)
    mc.stats<-numeric(n.mc)
    for(i in 1:n.mc){
      for(j in 1:n){
        mc.perm[j,]<-sample(possible.ranks[j,])
      }
      mc.stats[i]<-R.star.all(mc.perm)
    }
    for(i in 1:outp$num.comp){
      outp$p.val[i]<-mean(mc.stats>=outp$obs.stat[i])      
    }
  }  
  if(outp$method=="Asymptotic"){
    for(i in 1:outp$num.comp){
      adj<-outp$obs.stat[i]*(outp$n*outp$k*(outp$k+1)/6)^(-1/2)
      outp$p.val[i]<-pMaxCorrNor(adj,outp$k-1,rho=0.5)
    }
    
  }
  class(outp)<-"NSM3Ch7MCp"
  outp
}
