qKolSmirnLSA<-function(alpha){
  Qfun<-function(s){
    k<-c(-100:100)
    q<-(-1)^k*exp(-2*k^2*s^2)
    round(1-sum(q),4)
  }
  qalpha<-0
  upper<-0
  test<-seq(0.3, 2.39, 0.001)
  for(i in 1:length(test)){
    qalpha[i]<-test[i]
    upper[i]<-Qfun(test[i])
  }
  qtable<-matrix(c(qalpha,upper),ncol=2)
  pos<-which.min(abs(qtable[,2]-alpha))
  return(qtable[pos,1])
}