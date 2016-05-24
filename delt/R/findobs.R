findobs<-function(x,rec,ordobs,leftend,coordi){
#Finds which observations belong to given rectangle
#
#x on n*d havaintomatriisi
#rec on d*2-vector, sis kuvauksen uudesta osiosta
#ordobs is lkm-vector, points to the rows of x, ordered with respect
#  to i:th coordinate of x
#leftend >=0, <=lkm, integer, which pointers in ordobs point to the 
#  observations on the left rectangle
#coordi in 1:d, according to this coordinate pointers are ordered 
#
#Returns leftend in 0:lkm, endpoint of those belonging to rec
#
#We can start from leftend and proceed to the right hand side,
#because it was already checked that the observations on the left  
#side do not belong to rec.
#
lkm<-length(ordobs)
#
if (leftend<lkm){ #otherwise all obs already are on the left
   lcount<-0
   for (i in (leftend+1):lkm){
       ind<-ordobs[i]
       obs<-x[ind,]
       if (belongs(obs,rec,coordi))  lcount<-lcount+1  
   }
   leftend<-leftend+lcount
}
return(leftend)
}
