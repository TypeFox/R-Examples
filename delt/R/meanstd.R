meanstd<-function(apu){
#
wv<-dim(apu)[1]
tuldim<-dim(apu)[2]
cv<-matrix(0,tuldim,1)
cvstd<-matrix(0,tuldim,1)
#
for (lo in 1:tuldim){
   for (ek in 1:wv){
      cv[lo]<-cv[lo]+apu[ek,lo]
   }
   cv[lo]<-cv[lo]/wv
   for (ek in 1:wv){
      cvstd[lo]<-cvstd[lo]+(apu[ek,lo]-cv[lo])^2
   }
   cvstd[lo]<-sqrt(cvstd[lo]/wv)
   #
   #cv[lo]<-mean(apu[,lo])              
   #cvstd[lo]<-sd(apu[,lo])
}
#
return(list(cv=cv,cvstd=cvstd))
}
