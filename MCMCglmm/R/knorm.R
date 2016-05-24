knorm<-function(V,k){
  if(k%%2!=0){stop("k must be even")}
    k<-k/2
    m<-dim(V)[1]
    kfactorial<-factorial(k)
    KPM<-KPPM(m,k)
    MM<-matrix(0, m^k, m^k)
    for(h in 1:floor(k/2+1)){
      ah<-((kfactorial)^2)/(factorial(k-2*(h-1))*(factorial(h-1)^2)*(2^(2*(h-1))))
      Ir<-k-2*(h-1)
      Is<-h-1
      Ksigma<-1
      if(Ir!=0){
        for(i in 1:Ir){
          Ksigma<-kronecker(Ksigma, V)
        }
      }
      if(Is!=0){
        for(j in 1:Is){
          Ksigma<-kronecker(Ksigma, c(V)%*%t(c(V)))
        }
      }
      MM<-MM+ah*Ksigma
  }
  to.tensor(array(KPM%*%MM%*%KPM, dim=rep(m,k*2)))
}

