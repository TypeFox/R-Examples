skew<-function(x,method="unbiased"){
m.i<-function(x,i){
  n<-length(as.vector(x))
  m<-(1/n)*sum((x-mean(x))^i)
  m}

n<-length(as.vector(x))
if(method=="unbiased"){
skew<-(n/((n-1)*(n-2)))*sum(((x-mean(x))/sd(x))^3)}

if(method=="moments"){
skew<-m.i(x,3)/(m.i(x,2)^(3/2))}
skew
}
############################################################################
kurt<-function(x,method="unbiased"){
m.i<-function(x,i){
  n<-length(as.vector(x))
  m<-(1/n)*sum((x-mean(x))^i)
  m}

n<-length(as.vector(x))
if(method=="unbiased"){
kurtosis<-((n*(n+1))/((n-1)*(n-2)*(n-3)))*sum(((x-mean(x))/sd(x))^4)-
((3*((n-1)^2))/((n-2)*(n-3)))}

if(method=="moments"){
kurtosis<-m.i(x,4)/(m.i(x,2)^2)}

if(method=="excess"){
kurtosis<-(m.i(x,4)/(m.i(x,2)^2))-3}
kurtosis
}
