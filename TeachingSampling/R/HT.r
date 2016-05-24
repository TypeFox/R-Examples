HT<-function(y,Pik){
y<-t(as.matrix(y))
pik<-as.matrix(Pik)
HT<-y%*%(1/Pik)
HT
}