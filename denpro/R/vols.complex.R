vols.complex<-function(complex,dendat,meto="voltriangle")
{
# complex is lkm*(d+1) matrix
# dendat is n*d matrix

lkm<-dim(complex)[1]
vols<-matrix(0,lkm,1)

if (meto=="voltriangle"){
for (i in 1:lkm){
    ind<-complex[i,]
    simp<-dendat[ind,]
    vols[i]<-voltriangle(simp)
}}
else{
for (i in 1:lkm){
    ind<-complex[i,]
    simp<-dendat[ind,]
    vols[i]<-volsimplex(simp)
}}

return(vols)
}


