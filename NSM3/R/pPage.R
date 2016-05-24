pPage<-function(x,b=NA,trt=NA,method=NA, n.mc=10000){
  outp<-list()
  outp$stat.name<-"Page L"
  outp$n.mc<-n.mc  
  
  ties<-!length(unique(as.numeric(x)))==length(x)
  
  #If given a list, try to convert to a matrix. Each item 
  #in the list represents a column in the matrix.
  if(is.list(x)){x<-matrix(as.numeric(unlist(x)),ncol=length(x),byrow=F)}
  
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
    if(factorial(outp$k)^outp$n<=10000){
      method<-"Exact"
    }
    if(factorial(outp$k)^outp$n>10000){
      method<-"Monte Carlo"
    }
  }
  #####################################################################
  
  outp$method<-method
  
  L.calc<-function(x){
    return(sum((1:outp$k)*colSums(t(apply(x,1,rank)))))    
  }
  
  outp$obs.stat<-L.calc(x)
  
  possible.ranks<-t(apply(x,1,function(x) as.numeric(rank(x))))
  
  if(outp$method=="Exact"){
    possible.perm<-multCh7(possible.ranks)
    exact.dist<-numeric(factorial(outp$k)^outp$n)
    for(i in 1:factorial(outp$k)^outp$n){
      exact.dist[i]<-L.calc(possible.perm[,,i])
    }
  outp$p.val<-mean(exact.dist>=outp$obs.stat)      
  }

if(outp$method=="Monte Carlo"){
  mc.perm<-matrix(ncol=outp$k,nrow=outp$n)
  mc.stats<-numeric(n.mc)
  for(i in 1:n.mc){
    for(j in 1:n){
      mc.perm[j,]<-sample(possible.ranks[j,])
    }
    mc.stats[i]<-L.calc(mc.perm)
  }
  
  mc.vals<-sort(unique(mc.stats))
  mc.dist<-as.numeric(table(mc.stats))/n.mc
  outp$p.val<-mean(mc.dist>=outp$obs.stat) 
}

if(outp$method=="Asymptotic"){
    outp$stat.name<-"Page L*"
    outp$obs.stat<-(outp$obs.stat-((outp$n*outp$k*(outp$k+1)^2)/4))/(sqrt(outp$n*outp$k^2*(outp$k+1)*(outp$k^2-1)/144))
  } 
  outp$p.val<-1-pnorm(outp$obs.stat)

class(outp)<-"NSM3Ch7p"
outp
}