cfa.CVPRealz <-
function(sym.data,TFilas,TColumnas,TT,z) {
  ## z = eigenvector of matrix Z
  ## TT = Total
  NFil<-sym.data$N
  MCol<-sym.data$M
  aMin<-min(NFil,MCol)
  TablaDatos<-matrix(0,NFil,MCol) # To centers
  VPRealz<-matrix(0,aMin,aMin) # To Z
  if(MCol<=NFil) {   
    for(i in 1:aMin) {
       for(j in 1:aMin) { 
          VPRealz[i,j]<-abs(TT)*(sqrt(abs(TT))*sqrt(TColumnas[i]))*z[i,j]
       }
    }
  }
  else {	
     for(i in 1:aMin) {
       for(j in 1:aMin) { 
         VPRealz[i,j]<-abs(TT)*(sqrt(abs(TT))*sqrt(TFilas[i]))*z[i,j]
       }
     }
  }
  return(VPRealz)
}
