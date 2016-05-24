omaord2<-function(a,b){
#Jarjestaa vektorin a vektorin b mukaiseen jarjestykseen
#
#a and b are lnum-vectors
#
lnum<-length(a)  
orda<-a               #tahan oikea jarjestys
ordb<-b
i<-1 
while (i<=lnum){
   pienin<-omaind(b)
   ordb[i]<-b[pienin]
   orda[i]<-a[pienin]
   b[pienin]<-NA         #NA on plus aareton
   i<-i+1
}
return(orda)
}





