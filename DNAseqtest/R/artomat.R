artomat <-function(fobs){
n1<-fobs
nn=length(dim(fobs))
ee=NULL
for(i in 1:nn){
ee[[i]]=1:4
}
dimnames(n1)<-ee
w<-as.data.frame.table(n1)
w<-as.matrix(w)
w1<-w[,1:(dim(w)[2])-1]
w1<-matrix((as.numeric(w1)),dim(w)[1],dim(w)[2]-1)
w2<-cbind(as.numeric(w[,(dim(w)[2])]))
seq<-NULL
seq1<-NULL
for(i in 1:dim(w)[1]){
seq1<-rep(w1[i,],w2[i,])
seq<-c(seq,seq1)
}
fseq<-(matrix(seq,ncol=(dim(w1)[2]),byrow=T))
fseq
}
