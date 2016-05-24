tauK<-function(X){
 n<-sum(X)
 pi<-apply(X,1,sum)/n
 pj<-apply(X,2,sum)/n
 pij<-X/n
 s1<-sweep(pij^2,1,pi,"/")
 tau<-(sum(s1)-sum(pj^2))/(1-sum(pj^2))
 tau}