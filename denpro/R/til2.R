til2<-function(runi,curkosk,currecs,parimat,kosk){
#
blokki<-100
bloknum<-1
curpit<-blokki  # curpit<-lkm*curlkm #pahimillaan jokainen leikkaa jokaista
#
#if (dim(t(parimat))[1]==1) lkm<-1 else lkm<-length(parimat[1,]) 
lkm<-length(parimat[1,])    #kaiteiden maara
curlkm<-length(curkosk[,1]) #curlkm on edellisten kosketusten maara 
edkosk<-length(curkosk[1,]) #aikaisemmin haettiin kaikki leikkaukset
                            #edkosk:n kaiteen valilla  
d<-length(currecs[1,])/2
uuskosk<-matrix(0,curpit,kosk)     #matrix(0,lkm*curlkm,kosk) 
uusrecs<-matrix(0,curpit,2*d)      #matrix(0,lkm*curlkm,2*d)
#
ind<-0
for (i in 1:curlkm){  #kaydaan lapi kaikki nykyiset suorakaiteet
  vipu<-curkosk[i,edkosk] #uuden leikkaavan pitaa leikata esim viim curkosk:ssa
                  #to fasten the algorithm we do not consider every rec:
	          #only those who intersec with 1. in curkosk 
  ehdind<-1      
  ehdokas<-parimat[vipu,ehdind]  #ehdokkaat ovat ne jotka leikkaavat vipua
  while ((!is.na(ehdokas)) && (ehdind<=lkm)){ 
                  #kayd lapi ne jotka leikk vipua
                  #ehdokkaan pitaa leikata kaikkia muitakin 
                  #curkosk:n i:nnella rivilla olevia
    if (ehdokas>vipu){  #hetaan vain suuremmista kuin vipu   
      j<-1     
      touch<-TRUE
      olimuita<-FALSE
      while ((j<=(edkosk-1)) && (touch)){ #kayd lapi muut kuin vipu
        muu<-curkosk[i,j]
        if (!(ehdokas==muu)){
          olimuita<-TRUE
          curkoske<-parimat[ehdokas,]   #ne joihin ehdokas koskettaa
          touch<-onko(curkoske,muu)     #onko muu rivilla "curkoske"  
           #if (parimat2[ehdokas,muu]==0) touch<-FALSE
        }
              #jos ehdokas ja muu eivat kosketa ja ovat eri
	      #jos ehdokas=muu, niin parimat2[ehdokas,muu]=0  
        j<-j+1
      }
      if ((touch) && !(olimuita)) touch<-FALSE
      if (touch){  #jos ehdokas kosketti kaikkia muita
          ind<-ind+1   #lisataan uusien leikkausten laskuria
          if (ind<=curpit){  #jos ei tarvita uutta blokkia     
            uuskosk[ind,]<-c(curkosk[i,],ehdokas)  #???????
            uusrecs[ind,]<-leikkaa(currecs[i,],runi[ehdokas,])
          }
          else{
            bloknum<-bloknum+1
            uuspit<-bloknum*blokki
            apukosk<-matrix(0,uuspit,kosk)
            apurecs<-matrix(0,uuspit,2*d)
            apukosk[1:curpit,]<-uuskosk[1:curpit,]
            apurecs[1:curpit,]<-uusrecs[1:curpit,]
            apukosk[ind,]<-c(curkosk[i,],ehdokas)  #???????
            apurecs[ind,]<-leikkaa(currecs[i,],runi[ehdokas,])
            uuskosk<-apukosk
            uusrecs<-apurecs
            curpit<-uuspit
          }             
      }
    }  
    ehdind<-ehdind+1
    if (ehdind<=lkm) ehdokas<-parimat[vipu,ehdind] 
                                     #otetaan uusi vipua leikkaava
  }
}
if (ind>0){ 
       curkosk<-uuskosk[1:ind,]
       currecs<-uusrecs[1:ind,]
}     
return(list(ind=ind,currecs=currecs,curkosk=curkosk))
}








