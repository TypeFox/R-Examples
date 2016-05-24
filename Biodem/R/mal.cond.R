mal.cond <- function(PHI,N){
  k<-N/sum(N) ## the relative population of each k populaion on the total population of the area in study
  rmu<-PHI%*%k ## k is a list, coerced to vertical vector. Here I calculate the row wheight phi mean
  mu<-k%*%rmu ## k is now coerced to a linear vector. Here I calculated the overall mean phi
  az<-matrix(rep(rmu,length(rmu)),ncol=length(rmu))
  ax<-az+t(az)
  mu<-as.numeric(mu)
  r.mat<-(PHI+mu-ax)/(1-mu)
}
