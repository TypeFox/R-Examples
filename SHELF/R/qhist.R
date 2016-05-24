
qhist<-function(q,z,pz){
  approx(pz,z,q)$y
}