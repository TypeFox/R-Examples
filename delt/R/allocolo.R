allocolo<-function(mlkm,pg,mlkmpre,pgpre,colopre,coloind,colofall,paletti)
{

lkmpre<-mlkmpre$lkm   #previous number of modes
lkm<-mlkm$lkm         #current number of modes
modecolo<-matrix("",lkm,1)

# calculate distances

if (!is.null(colopre)){
dist<-matrix(NA,lkm,lkmpre)   #NA is infty
for (i in 1:lkm){
   for (j in 1:lkmpre){
      cent<-pg$center[,mlkm$modloc[i]]
      centpre<-pgpre$center[,mlkmpre$modloc[j]]
      dist[i,j]<-etais(cent[1:2],centpre[1:2])
   }
}
}

# allocate colors

if (is.null(colopre)){
for (k in 1:lkm){
  modecolo[k]<-paletti[k]
}
coloind<-lkm
}
else if (lkm>=lkmpre){
for (k in 1:lkmpre){
  minimi<-min(dist,na.rm=TRUE)
  argmin<-which(minimi==dist)[1]
  yind<-ceiling(argmin/lkm)      
  xind<-argmin-(yind-1)*lkm       #index for first color
  dist[xind,]<-NA
  modecolo[k]<-colopre[yind]
}
k<-lkmpre+1
while (k<=lkm){
  modecolo[k]<-paletti[coloind+k-lkmpre]
  k<-k+1
}
coloind<-lkm
}
else{                 #lkm<lkmpre
for (k in 1:lkm){
   minimi<-min(dist,na.rm=TRUE)
   argmin<-which(minimi==dist)[1]
   #cur<-omaindmat(dist)
   yind<-ceiling(argmin/lkm)      
   xind<-argmin-(yind-1)*lkm       #index for first color
   dist[xind,]<-NA
   modecolo[k]<-colopre[yind]
}
colofall<-lkm
}

return(list(modecolo=modecolo,coloind=coloind,colofall=colofall))
}




