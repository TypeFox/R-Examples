moodilkm<-function(vanhat)
{
#Lasketaan moodien lukumaara tiheyspuusta.
#Tiheyspuusta kaytettavissa vektori vanhat.
#Mikali solmu ei ole minkaan solmun vanhempi, se on lehti.

pit<-length(vanhat)
leima<-matrix(0,pit,1)
i<-1
while (i<=pit){
  solmu<-vanhat[i]
  leima[solmu]<-1
  i<-i+1
}
eimoodi<-sum(leima)
lkm<-(pit-eimoodi)
ykk<-rep(1,pit)
modnodes<-ykk-leima
#
moodiloc<-matrix(0,lkm,1)
ind<-1
for (i in 1:pit){
  if (modnodes[i]==1){
       moodiloc[ind]<-i
       ind<-ind+1
  }
}
#
return(list(lkm=lkm,modnodes=t(modnodes),modloc=t(moodiloc)))
}


