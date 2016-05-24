decom<-function(lsets,levels,links,alkublokki,blokki){
#Makes density tree,  edellyttaa jarjestyksen ???????
#
#lsets is levnum*atomnum-matrix
#levels is levnum-vector
#links is  atomnum*maxtouchnum-matrix, pointers to atoms,
#  for each atom we indicate which atoms it touches
#
#lsets matriisiin lisataan riveja, silla aikaisemmat rivit jakautuvat
#erillisiin osiinsa, lisataan levels:n vastaavat alkiot 
#
#Returns list(lsets,levels,parents)
# parents and levels are newlevnum-vectors
# lsets is newlevnum*atomnum
#
if (dim(t(lsets))[1]==1) levnum<-1 else levnum<-length(lsets[,1])
                                              #rows of lsets 
if (levnum==1) atomnum<-length(lsets) else
atomnum<-length(lsets[1,])     #maxalkio on n
newlevnum<-levnum*atomnum       #karkea arvio, jokainen tasojoukko 
#                               #voi periaatteessa jakautua n:aan osaan
newlsets<-matrix(0,alkublokki,atomnum) 
newlevels<-matrix(0,alkublokki,1)
parents<-matrix(0,alkublokki,1)
#
curblokki<-alkublokki
#
pino<-matrix(0,newlevnum,2) #1.col indeksi newlsets:iin, 2.col ind tasoon
a<-1
b<-2
#
if (levnum==1) levset<-lsets else
levset<-lsets[1,]                   #alin tasojoukko
rindeksit<-change(levset)            #change the representation
kumu<-declev(rindeksit,links,atomnum)  #jaetaan alin tasoj. osiin
if (dim(t(kumu))[1]==1) koko<-1 else koko<-length(kumu[,1]) #osien lkm
# Talletetaan osat 
newlsets[1:koko,]<-kumu   
newlevels[1:koko]<-levels[1]     #arvo toistuu 
efek<-koko                       #kirjataan uusien osien lkm
# Laitetaan kaikki osat pinoon
pino[1:koko,a]<-seq(1,koko)      #1,2,...,koko
pino[1:koko,b]<-1   #kaikki osat kuuluvat alimpaan tasojoukkoon
pinind<-koko        #indeksi pinoon
# 
if (levnum>1){  while (pinind>=1){
  ind<-pino[pinind,a]         #indeksi tasoon
  levind<-pino[pinind,b]      #ko tason korkeus
  pinind<-pinind-1            #otettu pinosta
  partlevset<-newlsets[ind,]
  higlevset<-lsets[levind+1,] #huom levind<levnum
  levset<-partlevset*higlevset
  if (sum(partlevset-levset)==0){  
          #jos leikkaus ei muuta, niin tasoj sailyy samana  
      newlevels[ind]<-levels[levind+1]  #poistetaan alempi osa      
      if (levind+1<levnum){ #jos ei olla korkeimmalla tasolla,laita pinoon
        pino[pinind+1,a]<-ind     
        pino[pinind+1,b]<-levind+1 #tasojouk taso on levind+1 
        pinind<-pinind+1       
      }
  }
  else if (sum(levset)>0){   #leikkaus ei tyhja
      rindeksit<-change(levset)
      kumu<-declev(rindeksit,links,atomnum) #jaet tasoj. osa osiin
      if (dim(t(kumu))[1]==1) koko<-1 else koko<-length(kumu[,1]) #osien lkm
      if ((efek+koko)>curblokki){   
        newlsets<-blokitus(newlsets,blokki)
        newlevels<-blokitus(newlevels,blokki)
        parents<-blokitus(parents,blokki)
        curblokki<-curblokki+blokki    
      }
      newlsets[(efek+1):(efek+koko),]<-kumu  #paivitetaan kumu tulokseen
      newlevels[(efek+1):(efek+koko),]<-levels[levind+1]   #arvo toistuu  
      parents[(efek+1):(efek+koko),]<-ind
      efek<-efek+koko
      if (levind+1<levnum){ #jos ei olla korkeimmalla tasolla,laita pinoon
        pino[(pinind+1):(pinind+koko),a]<-seq(efek-koko+1,efek)
        pino[(pinind+1):(pinind+koko),b]<-levind+1 #tasojouk taso on levind+1 
        pinind<-pinind+koko       
      }
  }
}} 
newlevels<-t(newlevels[1:efek])
newlsets<-newlsets[1:efek,]
parents<-t(parents[1:efek])
return(list(lsets=newlsets,levels=newlevels,parents=parents))
}








