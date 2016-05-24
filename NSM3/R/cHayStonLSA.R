cHayStonLSA<-function(alpha,k,delta=.001){
  
  ###Uses Hayter-Liu 1996
  our.grid<-seq(-8,8,delta)
  len<-length(our.grid)
  
  cutoffs<-rev(our.grid)
  for(iter in 1:len){
    init.grid<-pnorm(our.grid+cutoffs[iter])
  
    if(k>2){
      for(i in 3:k){
        new.grid<-cumsum(init.grid*dnorm(our.grid)*delta)+init.grid*(pnorm(our.grid+cutoffs[iter])-pnorm(our.grid))
        init.grid<-new.grid
      }
    }
  
    if(sum(dnorm(our.grid)*init.grid*delta)<=(1-alpha)){return(cutoffs[iter])}
  }
}