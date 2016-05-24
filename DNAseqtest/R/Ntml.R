Ntml <-
function(N,Ft){
s1<-length(dim(Ft))
s2<-length(Ft)
x<-array(0,c(rep(4,s1)))
x1<-0
ft1<-0
x[1]<-rbinom(1,N,Ft[1])
for(i in 2:(s2-1)){
x1<-x1+x[i-1]
ft1<-ft1+Ft[i-1]
x[i]<-rbinom(1,(N-x1),(Ft[i])/(1-ft1))
}
x1<-x1+x[s2-1]
ft1<-ft1+Ft[s2-1]
x[s2]<-N-x1
x
}
