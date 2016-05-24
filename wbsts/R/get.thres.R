get.thres <-
function(n, q=.95, r=100, scales=NULL){
  J<-round(log(n, 2))
  if(is.null(scales)) scales<-J-1:floor(J/2)
  max.scale<-length(scales)
  M<-NULL
  for(i in 1:r){
    x<-rnorm(n)
    ews<-ews.trans(x,scales=scales)
    m<-NULL
    for(l in 1:max.scale){
      z<-ews[,l]
      dis<-c(round((n-n/2^(scales[l])+1)):n)
      z<-z[-dis]; nz<-length(z)
      m<-c(m, max(abs(finner.prod.iter(z)))/1/log(n))
    }
    M<-rbind(M, m)
  }
  Sigma<-matrix(0, n, n)
  for(i in 1:n){ for(j in 1:n){
    Sigma[i,j]<-.3^abs(i-j)
  }}
  X<-mvtnorm:: rmvnorm(r, sigma=Sigma,method="chol")
  for(i in 1:r){
    x<-X[i,]
    ews<-ews.trans(x,scales=scales)
    m<-NULL
    for(l in 1:max.scale){
      z<-ews[,l]
      dis<-c(round((n-n/2^(scales[l])+1)):n)
      z<-z[-dis]; nz<-length(z)
      m<-c(m, max(abs(finner.prod.iter(z)))/1/log(n))
    }
    M<-rbind(M, m)
  }
  for(i in 1:n){ for(j in 1:n){
    Sigma[i,j]<-.6^abs(i-j)
  }}
  X<-mvtnorm:: rmvnorm(r, sigma=Sigma,method="chol")
  for(i in 1:r){
    x<-X[i,]
    ews<-ews.trans(x,scales=scales)
    m<-NULL
    for(l in 1:max.scale){
      z<-ews[,l]
      dis<-c(round((n-n/2^(scales[l])+1)):n)
      z<-z[-dis]; nz<-length(z)
      m<-c(m, max(abs(finner.prod.iter(z)))/1/log(n))
    }
    M<-rbind(M, m)
  }
  
  for(i in 1:n){ for(j in 1:n){
    Sigma[i,j]<-.9^abs(i-j)
  }}
  X<-mvtnorm:: rmvnorm(r, sigma=Sigma,method="chol")
  for(i in 1:r){
    x<-X[i,]
    ews<-ews.trans(x,scales=scales)
    m<-NULL
    for(l in 1:max.scale){
      z<-ews[,l]
      dis<-c(round((n-n/2^(scales[l])+1)):n) 
      z<-z[-dis]; nz<-length(z)
      m<-c(m, max(abs(finner.prod.iter(z)))/1/log(n))
    }
    M<-rbind(M, m)
  }
  return(apply(M, 2, function(x){quantile(x, q)}))
}
