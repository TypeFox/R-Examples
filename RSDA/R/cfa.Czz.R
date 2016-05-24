cfa.Czz <-
function(sym.data,TFilas,TColumnas,VPRealz,d) {
  NFil<-sym.data$N
  MCol<-sym.data$M
  aMin<-min(NFil,MCol)
  aMax<-max(NFil,MCol)
  X<-matrix(0,NFil,MCol) # To centers
  zz<-matrix(0,NFil,MCol) # To z
  for(i in 1:NFil) {
    for(j in 1:MCol) {
      X[i,j]<-(sym.var(sym.data,j)$var.data.vector[i,1]+
               sym.var(sym.data,j)$var.data.vector[i,2])/2
    }
  }  
  suma<-0
  if(MCol<=NFil) {   
    for(i in 1:aMin) {
      for(l in 1:aMax) {
        suma<-0
        for(s in 1:aMin) {
           suma<-suma+(X[l,s]/TColumnas[s])*VPRealz[s,i]
        }
        suma<-(1/sqrt(d[i]))*suma
        zz[l,i]<-suma;
      }
    }
  }
  else {
    for(i in 1:aMin) {
      for(l in 1:aMax) {
        suma<-0
        for(s in 1:aMin) {
          suma<-suma+(X[s,l]/TFilas[s])*VPRealz[s,i]
        }
        suma<-(1/sqrt(d[i]))*suma
        zz[l,i]=suma
      }
    }
  }
  return(zz)
}
