til1<-function(runi){
#
if (dim(t(runi))[1]==1) lkm<-1 else lkm<-length(runi[,1])
d<-length(runi[1,])/2
#
if (lkm>=2){
 parimat<-matrix(NA,lkm,lkm) #rivilla i lueteltu ne jotka leikk suorakaid i
                      # blokit tahan !!!!!!!!!!!!!!!!
 touchlkm<-matrix(0,lkm,1)   #kuinka monta kosketusta riv. i olevalle kaiteelle
# parimat2<-matrix(0,lkm,lkm) #rivilla i saralla j on 1 jos i ja j suork.leikk.
 l<-choose(lkm,2)
 curkosk<-matrix(0,l,2)      # blokit tahan !!!!!!!!!!!!!!!!
 currecs<-matrix(0,l,2*d)
 ind<-0
 for (i in 1:lkm){
   viite<-1 
   j<-i+1
   while (j<=lkm){
     ise<-leikkaa(runi[i,],runi[j,])
     if (!is.na(ise)){
       ind<-ind+1
       curkosk[ind,]<-c(i,j)
       currecs[ind,]<-ise 
       touchlkm[i]<-touchlkm[i]+1
       touchlkm[j]<-touchlkm[j]+1    
       parimat[i,touchlkm[i]]<-j
       parimat[j,touchlkm[j]]<-i
#       parimat2[i,j]<-1
#       parimat2[j,i]<-1
     }
     j<-j+1
   }
 }
}
if (ind==1){             #jos oli vain yksi leikkaus
  curkosk<-t(curkosk[1:ind,])
  currecs<-t(currecs[1:ind,])
  }
else if (ind>=2){      #jos oli useampi kuin yksi leikkaus
    curkosk<-curkosk[1:ind,]
    currecs<-currecs[1:ind,]
}
# supistetaan parimat
maxkosk<-max(touchlkm)
parimat<-parimat[,1:maxkosk]
if (maxkosk==1) parimat<-t(t(parimat))
return(list(ind=ind,curkosk=curkosk,currecs=currecs,parimat=parimat))
}






