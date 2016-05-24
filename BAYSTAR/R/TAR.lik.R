TAR.lik<-function(ay,p1,p2,ph.1,ph.2,sig.1,sig.2,lagd,thres,lagp1,lagp2,constant=1,thresVar){
n<-length(ay)
p<-max(max(lagp1),max(lagp2))+constant
p.1<-matrix(ph.1,nrow=p1+constant,ncol=1)   ## Build a matrix p.1 for ph.1
p.2<-matrix(ph.2,nrow=p2+constant,ncol=1)   ## Build a matrix p.2 for ph.2

if (!missing(thresVar)){
  if (length(thresVar) > n ){
  zt <- thresVar[1:n]
  cat("Using only first", n, "elements of threshold Variable\n")
  }
  else {zt<-thresVar}
lag.y<- zt[(p+1-lagd):(n-lagd)] 
}
else lag.y<-ay[(p+1-lagd):(n-lagd)]
yt<-ay[(p+1):n]
n1<-sum(lag.y<=thres); n2<-sum(lag.y>thres)  ## Count no. of observations for each regime

y.1<-matrix(yt[lag.y<=thres],ncol=1)     ## Arrange vector y and matrix X for each regime
y.2<-matrix(yt[lag.y>thres],ncol=1)

x.1<-matrix(NA,nrow=p1,ncol=n-p)      ## Arrange matrix X of regime 1
for (i in 1:p1){
x.1[i,]<-ay[(p-lagp1[i]+1):(n-lagp1[i])]}
if(p1 > 1) {
if (constant==1){
tx.1<-cbind(1,t(x.1[,lag.y<=thres]))}
else {tx.1<-t(x.1[,lag.y<=thres])}
}
if(p1 == 1) {
if (constant==1){
tx.1<-cbind(1,t(t(x.1[,lag.y<=thres])))}
else {tx.1<-t(t(x.1[,lag.y<=thres]))}
}
x.2<-matrix(NA,nrow=p2,ncol=n-p)      ## Arrange matrix X of regime 2
for ( i in 1:p2){
x.2[i,]<-ay[(p-lagp2[i]+1):(n-lagp2[i])]}

if(p2 > 1){
if (constant==1){
tx.2<-cbind(1,t(x.2[,lag.y>thres]))}
else {tx.2<-t(x.2[,lag.y>thres])}
}
if(p2 == 1){
if (constant==1){
tx.2<-cbind(1,t(t(x.2[,lag.y>thres])))}
else {tx.2<-t(t(x.2[,lag.y>thres]))}
}
ln.li<--((t(y.1-tx.1%*%p.1)%*%(y.1-tx.1%*%p.1))/(2*sig.1))-   ## The model log-likelihood function is
        ((t(y.2-tx.2%*%p.2)%*%(y.2-tx.2%*%p.2))/(2*sig.2))-   ## the sum of two normal distributions
        ((n1/2)*log(sig.1))-((n2/2)*log(sig.2))

return(ln.li)
}

