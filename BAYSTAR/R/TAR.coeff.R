TAR.coeff<-function(reg,ay,p1,p2,sig,lagd,thres,mu0,v0,lagp1,lagp2,constant=1,thresVar){
p<-max(max(lagp1),max(lagp2))+constant
n<- length(ay)

if (!missing(thresVar)){
  if (length(thresVar) > n ){
  zt <- thresVar[1:n]
  cat("Using only first", n, "elements of threshold Variable\n")
  }
  else zt<-thresVar
lag.y<- zt[(p+1-lagd):(n-lagd)] 
}
else lag.y<- ay[(p+1-lagd):(n-lagd)]

yt<- ay[(p+1):n]

if (reg==1){
ph<-rep(0.01,p1)
y.1<-matrix(yt[lag.y<=thres],ncol=1)  ## Arrange vector y of regime 1
x.1<-matrix(NA,nrow=p1,ncol=n-p)      ## Arrange matrix X of regime 1
for (i in 1:p1){
x.1[i,]<-ay[(p-lagp1[i]+1):(n-lagp1[i])]}

if(p1>1){
if (constant==1){
tx<-cbind(1,t(x.1[,lag.y<=thres]))
}
else {
tx<-t(x.1[,lag.y<=thres])
}
}
if(p1 == 1){
if (constant==1){
tx<-cbind(1,t(t(x.1[,lag.y<=thres])))
}
else {
tx<-t(t(x.1[,lag.y<=thres]))
}
}
yt<- matrix(yt[lag.y<=thres],ncol=1)


sigma<- (t(tx)%*%tx)/sig+v0    ## Variance of conditional posterior distribution
                               ## Mean of conditional posterior distribution
mu<- solve(sigma,((t(tx)%*%tx)/sig)%*%(solve((t(tx)%*%tx),t(tx)%*%yt))+v0%*%mu0)

ph<- rmvnorm(n = 1, mu, solve(sigma),method="chol")   ## Draw ph from the multivariate normal distribution

}

else {
ph<-rep(0.01,p2)
y.2<-matrix(yt[lag.y>thres],ncol=1)   ## Arrange vector y of regime 2
x.2<-matrix(NA,nrow=p2,ncol=n-p)      ## Arrange matrix X of regime 2
for ( i in 1:p2){
x.2[i,]<-ay[(p-lagp2[i]+1):(n-lagp2[i])]}

if(p2 > 1) {
if (constant==1){
tx<-cbind(1,t(x.2[,lag.y>thres]))
}
else {
tx<-t(x.2[,lag.y>thres])
}
}
if(p2 == 1){
if (constant==1){
tx<-cbind(1,t(t(x.2[,lag.y>thres])))
}
else {
tx<-t(t(x.2[,lag.y>thres]))
}
}
yt<- matrix(yt[lag.y>thres],ncol=1)

sigma<- (t(tx)%*%tx)/sig+v0     ## Variance of conditional posterior distribution
                                ## Mean of conditional posterior distribution
mu<- solve(sigma,((t(tx)%*%tx)/sig)%*%(solve((t(tx)%*%tx),t(tx)%*%yt))+v0%*%mu0)

ph<- rmvnorm(n=1,mu,solve(sigma),method="chol")   ## Draw ph from the multivariate normal distribution
}

return(ph)
}

