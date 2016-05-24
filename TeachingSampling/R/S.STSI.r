S.STSI<-function(S,Nh,nh)
{
S<-as.factor(S)
S<-as.factor(as.integer(S))
cum<-cumsum(nh)
sam<-matrix(0,sum(nh))
for(k in 1: length(nh)){
h<-which(S==k)
sam.h<-sample(Nh[k],nh[k])
if(k==1){
sam[1:nh[k]]<-h[sam.h]
}
if(k>1){
sam[(cum[k-1]+1):(cum[k])]<-h[sam.h]
}
}
sam
}

