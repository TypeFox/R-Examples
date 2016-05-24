digit<-function(luku,base){
#Gives representation of luku for system with base
#
#luku is a natural number >=0
#base is d-vector of integers >=2, d>=2, 
#base[d] tarvitaan vain tarkistamaan onko luku rajoissa
#
#Returns d-vector of integers.
#
#example: digit(52,c(10,10)), returns vector (2,5)
#
d<-length(base)
digi<-matrix(0,d,1)
jako<-matrix(0,d,1)
jako[d]<-base[1]
for (i in (d-1):1){
  jako[i]<-base[d-i+1]*jako[i+1]
}
vah<-0
for (i in 1:(d-1)){
  digi[i]<-floor((luku-vah)/jako[i+1]) #if digi[i]>base[i], then ERROR
  vah<-vah+digi[i]*jako[i+1]
}
digi[d]<-luku-vah
# annetaan vastaus kaanteisesti se 2354 annetaan c(4,5,3,2)
# talloin vastaavuus sailyy base:n kanssa 
#apu<-matrix(0,d,1)
#for (i in 1:d){
#  apu[i]<-digi[d-i+1]
#}
apu<-digi[d:1]
return(apu)
}
