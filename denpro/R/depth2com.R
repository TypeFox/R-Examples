depth2com<-function(dep,N){
#
d<-length(N)
logn<-log(N,base=2)
cusu<-cumsum(logn)
ind<-1
while ((ind<=d) && ((dep-cusu[ind])>0)){
  ind<-ind+1
}
direc<-min(ind,d)
if (direc==1){
  depind<-dep
}
else{
  depind<-dep-cusu[direc-1]
}
return(list(direc=direc,depind=depind))
}
