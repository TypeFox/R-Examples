til<-function(runi){
#
if (dim(t(runi))[1]==1) lkm<-1 else lkm<-length(runi[,1])
masses<-matrix(0,lkm,1)
#
masses[1]<-sum(massat(runi))    #kaiteiden massojen summa
#
if (lkm>=2){
 apu<-til1(runi)
 ind<-apu$ind
 curkosk<-apu$curkosk
 currecs<-apu$currecs
 parimat<-apu$parimat
 if (ind>0){ #jos oli parittaisia leikkauksia
    masses[2]<-sum(massat(currecs)) #parittaisten leikkausten massojen summa
    kosk<-3
    while (ind>1){
      write(ind,file="apu",append=TRUE)
      apu2<-til2(runi,curkosk,currecs,parimat,kosk)
      ind<-apu2$ind
      if (ind>0){
        currecs<-apu2$currecs
        curkosk<-apu2$curkosk
        masses[kosk]<-sum(massat(currecs))
      }
      kosk<-kosk+1
    }
 }
}
res<-0                 # res<-til3(masses)
for (i in 1:lkm){
  res<-res+(-1)^(i-1)*masses[i]
}                                 
return(res)
}
