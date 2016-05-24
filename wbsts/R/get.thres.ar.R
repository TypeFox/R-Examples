get.thres.ar <-
function(y, q=.95, r=100, scales=NULL){
  n=length(y)
  J<-round(log(n, 2))
  if(is.null(scales)) scales<-J-1:floor(J/2)
  max.scale<-length(scales)
  M<-NULL
  ar.est=ar(y)$ar
  X=matrix(0,r,n)
  for (i in 1:r)  X[i,] = arima.sim(n = n, list(ar = c(ar.est),sd = 1))
  for(i in 1:r){
    x<-X[i,]
    ews<-ews.trans(x,scales=scales)
    m<-NULL
    for(l in 1:max.scale){
      z<-ews[,l]
      dis<-c((n-n/2^l+1):n)  
      z<-z[-dis]; nz<-length(z)
      m<-c(m, max(abs(finner.prod.iter(z)))/1/(log(n)))
    }
    M<-rbind(M, m)
  }
  
  return(apply(M, 2, function(x){quantile(x, q)}))
}
