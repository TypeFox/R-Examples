pRangeNor<-function(x,k){
  r<-function(x,n){
    inner.int<-function(s){
      exp(-s^2)*(pnorm(s+x/2)-pnorm(s-x/2))^(n-2)
    }
    return(n*(n-1)*exp(-x^2/4)/(2*pi)*integrate(inner.int,-Inf,Inf)$value)
  }
  
  approx.dens<-test.grid<-seq(0,10,.001)
  test.grid<-round(test.grid,3)
  for(i in 1:length(test.grid)){
    approx.dens[i]<-r(test.grid[i],k)*.001
  }
  approx.dens=approx.dens/sum(approx.dens)
  
  upper.tails<-rev(cumsum(rev(approx.dens)))
  upper.tails[test.grid==round(x,3)]
}
