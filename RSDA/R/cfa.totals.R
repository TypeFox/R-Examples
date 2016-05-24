cfa.totals <-
function(sym.data) {
  Total<-0  
  TotalMin<-0
  TotalMax<-0	
  N<-sym.data$N
  M<-sym.data$M
  TotalFilas<-rep(0,N)
  TotalColumnas<-rep(0,M)
  TotalFilasMin<-rep(0,N)
  TotalColumnasMin<-rep(0,M)
  TotalFilasMax<-rep(0,N)
  TotalColumnasMax<-rep(0,M)
  A<-matrix(0,N,M) # To centers
  CMin<-matrix(0,N,M)
  CMax<-matrix(0,N,M)  
  for(i in 1:N) {
    for(j in 1:M) {
      CMin[i,j]<-sym.var(sym.data,j)$var.data.vector[i,1]
      CMax[i,j]<-sym.var(sym.data,j)$var.data.vector[i,2]                  
      A[i,j]<-(sym.var(sym.data,j)$var.data.vector[i,1]+
               sym.var(sym.data,j)$var.data.vector[i,2])/2
    }
  }
  # Centers Totals
  for(i in 1:N)  {
    TotalFilas[i]<-sum(A[i,])
  }
  for(j in 1:M)  {
    TotalColumnas[j]<-sum(A[,j])
  }
  Total<-sum(TotalFilas)
  # To minimus
  for(i in 1:N)  {
    TotalFilasMin[i]=sum(CMin[i,])
  }
  for(j in 1:M)  {
    TotalColumnasMin[j]<-sum(CMin[,j])
  }
  TotalMin<-sum(TotalFilasMin)
  # To maximuns
  for(i in 1:N)  {
    TotalFilasMax[i]=sum(CMax[i,])
  }
  for(j in 1:M)  {
    TotalColumnasMax[j]<-sum(CMax[,j])
  }
  TotalMax<-sum(TotalFilasMax)  
  return(list(Total=Total,TotalMin=TotalMin,TotalMax=TotalMax,TotalRows=TotalFilas,
              TotalColumns=TotalColumnas,TotalRowsMin=TotalFilasMin,
              TotalColumnsMin=TotalColumnasMin,TotalRowsMax=TotalFilasMax,
              TotalColumnsMax=TotalColumnasMax))
}
