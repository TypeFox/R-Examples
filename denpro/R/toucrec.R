toucrec<-function(atoms,alkublokki,blokki){
#Finds which atoms touch each other
#
#items is atomnum*(2*d)-matrix
#alkublokki is an estimate to the maximum number of touches
#
#Returns links, atomnum*maxtouches-matrix
#
if (dim(t(atoms))[1]==1) m<-1 else m<-length(atoms[,1]) #m is number of atoms
len<-alkublokki
links<-matrix(NA,m,len)
maara<-matrix(0,m,1)
# merkitaan kosketukset linkit-matriisiin
i<-1
while (i<=m){
  j<-i+1
  while (j<=m){
    rec1<-atoms[i,]
    rec2<-atoms[j,]
    crit<-touch(rec1,rec2)
    if (crit){ #jos suorakaiteet koskettavat
        maari<-maara[i]+1
        maarj<-maara[j]+1
        if ((maari>len) || (maarj>len)){
            links<-blokitus2(links,blokki)
            len<-len+blokki
        }
        links[i,maari]<-j
        maara[i]<-maari
        links[j,maarj]<-i
        maara[j]<-maarj         
    }
    j<-j+1 
  }
  i<-i+1
} 
return(links)
}






