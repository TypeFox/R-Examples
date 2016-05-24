makeorder<-function(obspoint,x,coordi){
#Orders obspoint according to coordi:th coordinate of x
#
#obspoint is lkm-vector: pointers to rows of x
#x is n*d-matrix
#coordi is in 1:d
#
lkm<-length(obspoint)
#
redu<-matrix(0,lkm,1)
for (i in 1:lkm){
  obsind<-obspoint[i]
  redu[i]<-x[obsind,coordi]
}
ordobs<-matrix(0,lkm,1)
for (i in 1:lkm){
   pienin<-omaind(redu)
   ordobs[i]<-obspoint[pienin]
   redu[pienin]<-NA  #NA on plus aareton
}
return(ordobs)
}

