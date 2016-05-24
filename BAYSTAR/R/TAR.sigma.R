TAR.sigma <-function(reg,ay,thres,lagd,p1,p2,ph,v,lambda,lagp1,lagp2,constant=1,thresVar){
n<-length(ay)             ## Total no. of observations
p<-max(max(lagp1),max(lagp2))+constant
yt<- ay[(p+1):n]

if (!missing(thresVar)){
  if (length(thresVar) > n ){
  zt <- thresVar[1:n]
  cat("Using only first", n, "elements of threshold Variable\n")
  }
  else zt<-thresVar
lag.y<- zt[(p+1-lagd):(n-lagd)] 
}
else lag.y<- ay[(p+1-lagd):(n-lagd)]


if (reg==1){

m<- sum(lag.y<=thres)
y<- matrix(yt[lag.y<=thres],ncol=1)
x.1<-matrix(NA,nrow=p1,ncol=n-p)      ## Arrange matrix X of regime 1
for (i in 1:p1){
x.1[i,]<-ay[(p-lagp1[i]+1):(n-lagp1[i])]}

if(p1 > 1){
if (constant==1){
tx<-cbind(1,t(x.1[,lag.y<=thres]))}
else {tx<-t(x.1[,lag.y<=thres])}
}
if(p1 == 1){
if (constant==1){
tx<-cbind(1,t(t(x.1[,lag.y<=thres])))}
else {tx<-t(t(x.1[,lag.y<=thres]))}
}

phi<- matrix(ph,nrow=p1+constant,ncol=1)

s2<- (t(y-tx%*%phi)%*%(y-tx%*%phi))/m   ## Compute sample variance
}

else{
m<- sum(lag.y>thres)
phi<- matrix(ph,nrow=p2+constant,ncol=1)
y<- matrix(yt[lag.y>thres],ncol=1)
x.2<-matrix(NA,nrow=p2,ncol=n-p)      ## Arrange matrix X of regime 2
for ( i in 1:p2){
x.2[i,]<-ay[(p-lagp2[i]+1):(n-lagp2[i])]}
if(p2 > 1){
if (constant==1){
tx<-cbind(1,t(x.2[,lag.y>thres]))}
else {tx<-t(x.2[,lag.y>thres])}
}
if(p2 == 1){
if (constant==1){
tx<-cbind(1,t(t(x.2[,lag.y>thres])))}
else {tx<-t(t(x.2[,lag.y>thres]))}
}
s2<- (t(y-tx%*%phi)%*%(y-tx%*%phi))/m   ## Compute sample variance
}

shape<- (v+m)/2              ## Set shape parameter of Gamma distribution
rate<- (v*lambda+m*s2)/2     ## Set rate parameter of Gamma distribution
sigma<- 1/rgamma(1, shape=shape, rate=rate)  ## Draw sigma^2 from an Inverse-Gamma distribution.

return(sigma)
}

