cfa.MatrixZ <-
function(sym.data,TFilas,TColumnas) {
  Fil<-sym.data$N
  Col<-sym.data$M
  a<-min(Fil,Col)
  TablaDatos<-matrix(0,Fil,Col) # To centers
  A<-matrix(0,a,a) # To Z
  for(i in 1:Fil) {
    for(j in 1:Col) {
      TablaDatos[i,j]<-(sym.var(sym.data,j)$var.data.vector[i,1]+
                 sym.var(sym.data,j)$var.data.vector[i,2])/2
    }
  }  
  if(a==Col) {  
     for(j in 1:a) {
       for(l in j:a) {
         suma<-0
         for(i in 1:Fil) { 
           suma<-suma+((TablaDatos[i,j])*(TablaDatos[i,l]))/(TFilas[i])
         }
         suma<-suma*(1/sqrt(abs((TColumnas[j])*(TColumnas[l]))))
         A[j,l]<-suma
         if(i!=j)
           A[l,j]<-suma
       }
    }
  }
  else {
    for(j in 1:a) {
      for(l in j:a) {
        suma<-0
        for(i in 1:Col) { 
          suma<-suma+((TablaDatos[j,i])*(TablaDatos[l,i]))/(TColumnas[i])
        }
        suma<-suma*(1/sqrt(abs((TFilas[j])*(TFilas[l]))))
        A[j,l]<-suma
        if(i!=j)
          A[l,j]<-suma
      }
    }
  }
  return(A)
}
