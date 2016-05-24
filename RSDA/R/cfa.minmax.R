cfa.minmax <-
function(sym.data,TFilas,TFilasMin,TFilasMax,TColumnas,
                   TColumnasMin,TColumnasMax,Total,VP,VPzz) {
  n<-sym.data$N
  m<-sym.data$M
  aMin<-min(n,m)
  X<-matrix(0,n,m) # To centers
  XMin<-matrix(0,n,m) 
  XMax<-matrix(0,n,m) 
  A<-matrix(0,n+m,aMin-1) 
  Min<-matrix(0,n+m,aMin-1) 
  Max<-matrix(0,n+m,aMin-1) 
  for(i in 1:n) {
    for(j in 1:m) {
      XMin[i,j]<-sym.var(sym.data,j)$var.data.vector[i,1]
      XMax[i,j]<-sym.var(sym.data,j)$var.data.vector[i,2]                  
      X[i,j]<-(sym.var(sym.data,j)$var.data.vector[i,1]+
                 sym.var(sym.data,j)$var.data.vector[i,2])/2
    }
  }
  for(j in 1:m) {
    for(k in 1:n) {
      X[k,j]<-X[k,j]/Total
      XMin[k,j]<-XMin[k,j]/Total
      XMax[k,j]<-XMax[k,j]/Total
    }
  }  
  if(m<=n) {   
    for(alfa in 2:aMin) { 
        for(i in 1:n) {
          suma<-0
          SumaA<-0
          SumaB<-0
          SumaC<-0
          SumaD<-0
          for(j in 1:m) {
            suma=suma+(X[i,j]*VP[j,alfa])/TColumnas[j]
            if(VP[j,alfa]<0) {
              SumaA<-SumaA+(XMax[i,j]*VP[j,alfa])/TColumnas[j]
              SumaC<-SumaC+(XMin[i,j]*VP[j,alfa])/TColumnas[j]
            }
            if(VP[j,alfa]>0) {
              SumaB<-SumaB+(XMin[i,j]*VP[j,alfa])/TColumnas[j]
              SumaD<-SumaD+(XMax[i,j]*VP[j,alfa])/TColumnas[j]
            }
          }
          suma<-suma/TFilas[i]
          SumaA<-SumaA/TFilas[i]
          SumaB<-SumaB/TFilas[i]
          SumaC<-SumaC/TFilas[i]
          SumaD<-SumaD/TFilas[i]
          A[i,(alfa-1)]<-suma
          Min[i,(alfa-1)]<-(SumaA+SumaB)
          Max[i,(alfa-1)]<-(SumaC+SumaD)
        }
    }
    for(alfa in 2:aMin) { 
      for(j in 1:m) {
        suma<-0
        SumaA<-0
        SumaB<-0
        SumaC<-0
        SumaD<-0
        for(i in 1:n) {          
          suma<-suma+(X[i,j]*VPzz[i,alfa])/TFilas[i]
          if(VPzz[i,alfa]<0) {
            SumaA<-SumaA+(XMax[i,j]*VPzz[i,alfa])/TFilas[i]
            SumaC<-SumaC+(XMin[i,j]*VPzz[i,alfa])/TFilas[i]
          }
          if(VPzz[i,alfa]>0) {
            SumaB<-SumaB+(XMin[i,j]*VPzz[i,alfa])/TFilas[i]
            SumaD<-SumaD+(XMax[i,j]*VPzz[i,alfa])/TFilas[i]
          }
        }
        suma<-suma/TColumnas[j]
        SumaA<-SumaA/TColumnas[j]
        SumaB<-SumaB/TColumnas[j]
        SumaC<-SumaC/TColumnas[j]
        SumaD<-SumaD/TColumnas[j]
        A[(n+j),(alfa-1)]<-suma
        Min[(n+j),(alfa-1)]<-(SumaA+SumaB)
        Max[(n+j),(alfa-1)]<-(SumaC+SumaD)
      }
    }
  }
  else {	
    for(alfa in 2:aMin) { 
      for(j in 1:m) {
        suma<-0
        SumaA<-0
        SumaB<-0
        SumaC<-0
        SumaD<-0
        for(i in 1:n) {
          suma=suma+(X[i,j]*VP[i,alfa])/TFilas[i]
          if(VP[i,alfa]<0) {
            SumaA<-SumaA+(XMax[i,j]*VP[i,alfa])/TFilas[i]
            SumaC<-SumaC+(XMin[i,j]*VP[i,alfa])/TFilas[i]
          }
          if(VP[i,alfa]>0) {
            SumaB<-SumaB+(XMin[i,j]*VP[i,alfa])/TFilas[i]
            SumaD<-SumaD+(XMax[i,j]*VP[i,alfa])/TFilas[i]
          }
        }
        suma<-suma/TColumnas[j]
        SumaA<-SumaA/TColumnas[j]
        SumaB<-SumaB/TColumnas[j]
        SumaC<-SumaC/TColumnas[j]
        SumaD<-SumaD/TColumnas[j]
        A[j,(alfa-1)]<-suma
        Min[j,(alfa-1)]<-(SumaA+SumaB)
        Max[j,(alfa-1)]<-(SumaC+SumaD)
      }
    }
    for(alfa in 2:aMin) { 
      for(i in 1:n) {
        suma<-0
        SumaA<-0
        SumaB<-0
        SumaC<-0
        SumaD<-0
        for(j in 1:m) {
          suma<-suma+(X[i,j]*VPzz[j,alfa])/TColumnas[j]
          if(VPzz[j,alfa]<0) {
            SumaA<-SumaA+(XMax[i,j]*VPzz[j,alfa])/TColumnas[j]
            SumaC<-SumaC+(XMin[i,j]*VPzz[j,alfa])/TColumnas[j]
          }
          if(VPzz[j,alfa]>0) {
            SumaB<-SumaB+(XMin[i,j]*VPzz[j,alfa])/TColumnas[j]
            SumaD<-SumaD+(XMax[i,j]*VPzz[j,alfa])/TColumnas[j]
          }
        }
        suma<-suma/TFilas[i]
        SumaA<-SumaA/TFilas[i]
        SumaB<-SumaB/TFilas[i]
        SumaC<-SumaC/TFilas[i]
        SumaD<-SumaD/TFilas[i]
        A[(m+i),(alfa-1)]=suma
        Min[(m+i),(alfa-1)]<-(SumaA+SumaB)
        Max[(m+i),(alfa-1)]<-(SumaC+SumaD)
      }
    }
  }
  return(list(Centers=A,Min=Min,Max=Max))
}
