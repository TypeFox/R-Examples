montecarlo.complex<-function(dendat,complex,rho,M,seed=1)
{
# dendat on n*d matriisi
n<-dim(dendat)[1]
d<-dim(dendat)[2]
lkm<-dim(complex)[1]

# create Monte Carlo sample
  ota<-c(complex)
  dendat2<-dendat[ota,]
  lows<-matrix(0,d,1)
  higs<-matrix(0,d,1)
  for (i in 1:d){
      ma<-max(dendat2[,i])
      mi<-min(dendat2[,i])
      lows[i]<-mi
      higs[i]<-ma
  }
  set.seed(seed)
  varia<-matrix(runif(M*d),M,d)
  for (i in 1:d) varia[,i]<-varia[,i]*(higs[i]-lows[i])+lows[i]

# laske kuinka monta joukossa

count<-matrix(0,M,1)
for (i in 1:M){
    point<-varia[i,]
    sisalla<-0
    j<-1
    while ( (j<=lkm) && (sisalla==0) ){
        simpl<-complex[j,]
        simple<-dendat[simpl,]
        sisalla2<-is.inside.simp.bary(point,simple)
        #sisalla2<-is.inside.simp.long(point,simple,rho)
        #sisalla2<-is.inside.simp(point,simple,rho)
        if (sisalla2==1){ 
               sisalla<-1
               count[i]<-1
        }
        j<-j+1
    }
}

# lakse tilavuus

  volurec<-1
  for (i in 1:d) volurec<-volurec*(higs[i]-lows[i])
  volu<-volurec*sum(count)/M

return(volu)
}

