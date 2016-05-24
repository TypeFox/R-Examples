

phist<-function(x,z,pz){
  z<-c(-exp(100), z, exp(100))
  pz<-c(0,pz,1)
  approx(z,pz,x)$y
}

