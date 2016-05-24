dentree<-function(treeseq,dendat,leafnum=NULL){
#Returns a tree from a sequence of trees
#
#if leafnum is NULL, we use empirical smoothing parameter selection
#
n<-length(dendat[,1])
#poistettavat<-treeseq$leafs
#if (dim(t(poistettavat))[1]==1) maxrem<-1
#else maxrem<-length(poistettavat[,1])      
#
if (is.null(leafnum)){   #empirical smoothing param. selection
  ri<-riskesti(treeseq,n)
  indeksi<-omaind(ri)
#  indeksit<-detsikko(treeseq$leafs,indeksi) #Haet haluttu puu puit jonosta
#  alipuu<-poistamon(treeseq$tree,indeksit)
}
else{
  indeksi<-detsi(treeseq$info[,1],leafnum)
#  indeksit<-detsikko(treeseq$leafs,indeksi)
  #indeksit<-t(t(treeseq$leafs))[indeksi,] 
#  alipuu<-poistamon(treeseq$tree,indeksit)  #alipuu<-dpoimi(puuseq,indeksi)
}
#
# Tehdaan binaaripuusta paloittain vakio
#
xlkm<-length(dendat[1,])
epsi<-0.01
#kanta<-kantaja(dendat,epsi)
#palvak<-teeositus(alipuu,xlkm,kanta)
#
# Tehdaan paloittain vakiosta tiheyspuu
#
#return(list(binpuu=alipuu,palvak=palvak))
alipuu<-NULL
return(alipuu)
}



