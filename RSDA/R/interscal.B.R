interscal.B <-
function(sym.data) {
  MDis<-interscal.delta(sym.data)
  Fil<-dim(MDis)[1] # Fil is 2m
  suma1<-0
  suma2<-0
  suma3<-0
  B<-matrix(0,Fil,Fil)
  for(r in 1:Fil) {
    for(s in 1:Fil) {
      suma3<-suma3+(MDis[r,s]*MDis[r,s])
    }
  }
  for(i in 1:Fil) {
    for(j in 1:Fil) {
      suma1<-0
      for(r in 1:Fil) {
        suma1=suma1+(MDis[r,j]*MDis[r,j])
      }
      suma2<-0
      for(s in 1:Fil) {
        suma2<-suma2+(MDis[i,s]*MDis[i,s])
      }
      B[i,j]<-(-1/2)*(MDis[i,j]*MDis[i,j]-(1/Fil)*suma1-(1/Fil)*suma2+((1/Fil)^2)*suma3)
    }
  }
  return(B)
}
