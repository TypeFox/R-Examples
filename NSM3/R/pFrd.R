pFrd<-function(x,b=NA,trt=NA,method=NA, n.mc=10000){
  outp<-list()
  outp$stat.name<-"Friedman, Kendall-Babington Smith S"
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
  
  S.calc<-function(x){
    sum.squares<-sum(colSums(t(apply(x,1,rank)))^2)
    return(round(12/(n*k*(k+1))*sum.squares-3*n*(k+1),10))  
  }
  
  outp$obs.stat<-S.calc(x)

  possible.ranks<-t(apply(x,1,function(x) as.numeric(rank(x))))
  
  if(outp$method=="Exact"){
    if(!ties){
      phi<-function(full){
      mat<-full[,-1]
      sort<-t(apply(mat,1,sort))
      uniq<-unique(sort,MARGIN=1)
      sort<-cbind(full[,1],sort)
      counts<-numeric(nrow(uniq))
      for(i in 1:length(counts)){
        counts[i]<-sum(apply(sort,1,function(x,y) x[1]*identical(x[-1],y), y=uniq[i,]))
      }
      return(as.matrix(cbind(counts[],uniq)))
    }
  
    update<-function(full,original.with.ranks){
      mat<-full[,-1]
      original<-original.with.ranks[,-1]
      output<-matrix(nrow=dim(original)[1]*max(nrow(full),1),ncol=min(dim(mat)[2],length(mat))+1)
      count<-1
      for(i in 1:max(dim(mat)[1],1)){
        if(max(dim(mat)[1],1)==1){
          for(j in 1:max(dim(original)[1],1)){
            output[count,]<-c(full[1],mat+original[j,])
            count<-count+1
          }
        }
        if(max(dim(mat)[1],1)!=1){
          for(j in 1:max(dim(original)[1],1)){
            output[count,]<-c(full[i,1],mat[i,]+original[j,])	
            count<-count+1
          }
        }
      }
      return(output)
    }

    exact.friedman.dist<-function(k,n,STATISTIC){
      initial<-cbind(rep(1,factorial(k)),multComb(rep(1,k)))
      if(nrow(initial)!=factorial(k)){
        return("Error!")
      }
      new<-update(phi(initial),initial)
      for(i in 1:(n-2)){
        new<-phi(update(new,initial))
      }
  
      sum.squares<-apply(new[,-1]^2,1,sum)
      statistic<-round(12/(n*k*(k+1))*sum.squares-3*n*(k+1),4)
      test.dist<-cbind(new[,1]/sum(new[,1]),statistic)
      p.value<-sum(test.dist[statistic>=STATISTIC,1])
      return(p.value)
    }
  
    outp$p.val<-exact.friedman.dist(outp$k,outp$n,outp$obs.stat)
  }
  
  if(ties){
    possible.perm<-multCh7(possible.ranks)
    exact.dist<-numeric(factorial(outp$k)^outp$n)
    for(i in 1:(factorial(outp$k)^outp$n)){
      exact.dist[i]<-S.calc(possible.perm[,,i])
    }
    outp$p.val<-mean(exact.dist>=outp$obs.stat)      
  }
  }

  if(outp$method=="Monte Carlo"){
    mc.perm<-matrix(ncol=outp$k,nrow=outp$n)
    mc.stats<-numeric(n.mc)
    for(i in 1:n.mc){
      for(j in 1:n){
        mc.perm[j,]<-sample(possible.ranks[j,])
      }
      mc.stats[i]<-S.calc(mc.perm)
    }
    
    mc.vals<-sort(unique(mc.stats))
    mc.dist<-as.numeric(table(mc.stats))/n.mc
    outp$p.val<-mean(mc.dist>=outp$obs.stat) 
  }
  
  if(outp$method=="Asymptotic"){
    if(ties){
      tie.groups<-as.numeric(unlist(apply(x,1,function(x) as.numeric(table(x)))))
      adj.size<-1/(outp$k-1)*(sum(tie.groups^3)-outp$n*outp$k)
      outp$stat.name<-"Friedman, Kendall-Babington Smith S'"
      outp$obs.stat<-outp$obs.stat*(outp$n*outp$k*(outp$k+1))/((outp$n*outp$k*(outp$k+1))-adj.size)
    } 
    outp$p.val<-1-pchisq(outp$obs.stat,outp$k-1)
  }
  
  class(outp)<-"NSM3Ch7p"
  outp
}
