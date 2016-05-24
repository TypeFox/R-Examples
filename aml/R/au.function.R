
.grlik2<- function(theta,residu, hd){
## this function is used by amltest() to calcuate the gradient of the likelihood function

   nn<-length(residu)
   diff1<-nn/theta[1]-residu%*%diag(1/(hd+theta[2]))%*%t(residu)/theta[1]^2
   diff2<-sum(1/(hd+theta[2]))-residu%*%diag(1/(hd+theta[2])^2)%*%t(residu)/theta[1]
   res<-c(diff1,diff2)
   return(res)
}

.lik2<-function(theta, residu, hd){
## this function is used by almtest() to evaludate the likelihood

  nn<-length(residu)
  res<-nn*log(theta[1])+sum(log(hd+theta[2]))+residu%*%diag(1/(hd+theta[2]))%*%t(residu)/theta[1]
  return(res)
}




.nchoose2<- function(n,x){
##  function used by amltest() in calculating EBIC
  haf<- floor(n/2)
  if (x<= haf) res<-choose(n,x)
  if (x> haf) res<- choose(n,haf)
  return(res)
}

.calseq<-function(x){
##  function used by amltest() in counting nonzero elements
  res<-sum(x!=0)
  return(res)
}



