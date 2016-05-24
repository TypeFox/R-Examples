alpha.complex<-function(complex,dendat,alpha)
{
M<-dim(complex)[1]
n<-dim(dendat)[1]
d<-dim(dendat)[2]  # d<-dim(complex)[2]-1  

acomplex<-matrix(0,M,d+1)
lkm<-0
for (m in 1:M){
    simindex<-complex[m,]
    simplex<-dendat[simindex,]

    tulos<-0
    i<-1
    while ((i<=d) && (tulos==0)){
       v1<-simplex[i,]
       j<-i+1
       while ((j<=(d+1)) && (tulos==0)){
         v2<-simplex[j,]
         etais2<-sum((v1-v2)^2)
         if (etais2>alpha^2) tulos<-1
         j<-j+1
       }
       i<-i+1
    }
    if (tulos==0){ 
       lkm<-lkm+1
       acomplex[lkm,]<-complex[m,]
    }
}
acomplex<-acomplex[1:lkm,]

return(acomplex)
}

