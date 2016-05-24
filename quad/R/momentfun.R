momentfun <-
function(Px,Py,n,mycoef){

f.1=mycoef[[5]]
f.2=mycoef[[6]]
f.3=mycoef[[7]]
f.4=mycoef[[8]]

v.1=mycoef[[9]]
v.2=mycoef[[10]]
v.3=mycoef[[11]]
v.4=mycoef[[12]]
 
 # Px and Py are lists of linear combinations

 ###########
 ### r=1 ###
 ###########

 w<-"r1.p1"; f<-f.1; v<-v.1

 num<-Px[[w]]*Py[[w]]*f
 denom<-prod(n:(n-v+1))

 first<-num/denom

 ###########
 ### r=2 ###
 ###########

 w<-c("r2.p1","r2.p2","r2.p3"); f<-f.2; v<-v.2

 num<-unlist(Px[w])*unlist(Py[w])*f
 denom<-rep(NA,length(v))
 for(i in 1:length(v)){
  denom[i]<-prod(n:(n-v[i]+1))
  }
 
 second<-sum(num/denom)
 
 ###########
 ### r=3 ###
 ###########

 w<-c("r3.p1","r3.p2","r3.p3","r3.p4","r3.p5","r3.p6","r3.p7","r3.p8"); f<-f.3; v<-v.3

 num<-unlist(Px[w])*unlist(Py[w])*f
 denom<-rep(NA,length(v))
 for(i in 1:length(v)){
  denom[i]<-prod(n:(n-v[i]+1))
  }
 
 third<-sum(num/denom)

 ###########
 ### r=4 ###
 ###########

 w<-rep("",length(f.4))
 for(i in 1:length(w)){ w[i]<-paste("r4.p",i,sep="") }; f<-f.4; v<-v.4

 num<-unlist(Px[w])*unlist(Py[w])*f
 denom<-rep(NA,length(v))
 for(i in 1:length(v)){
  denom[i]<-prod(n:(n-v[i]+1))
  }
 
 fourth<-sum(num/denom)

 return(list(first=first,second=second,third=third,fourth=fourth))
 }
