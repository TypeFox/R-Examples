allosimp<-function(vec,indvec,colot){
#
#minimi<-min(vec,na.rm=TRUE)
#cur<-which(vec==min(vec,na.rm=TRUE)) 
#
epsi<-0.0000001
reslen<-length(vec)
res<-matrix("",reslen,1)
#
len<-length(indvec)
for (i in 1:len){
   inde<-indvec[i]
   #vektori<-which(vec==vec[inde])
   #vektori<-which(((vec<=vec[inde]+epsi)&&(vec>=vec[inde]-epsi)))
   apu<-matrix(FALSE,reslen,1)
   for (l in 1:reslen){
     if ((vec[l]<=vec[inde]+epsi)&&(vec[l]>=vec[inde]-epsi)) apu[l]<-TRUE
   }
   vektori<-which(apu)
   len2<-length(vektori)
   for (j in 1:len2){
      sijo<-vektori[j]
      res[sijo]<-colot[i]
   }
}
#
return(res)
}
