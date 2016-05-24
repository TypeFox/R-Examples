interscal.delta <-
function(sym.data) {
  res<-interval.dist(sym.data,distance='interscal')
  CMin<-res$min.matrix
  CMax<-res$max.matrix  
  Fil<-sym.data$N
  D2M<-matrix(0,2*Fil,2*Fil)
  for(i in 1:Fil) {
    for(j in i:Fil) {
      if(i==j) {
        if(j!=Fil) {
          D2M[2*i-1,2*j-1]<-0
          D2M[2*i,2*j]<-0
          D2M[2*i-1,2*j]<-CMax[i,j]
          D2M[2*j,2*i-1]<-D2M[2*i-1,2*j]
        }
        else {
          D2M[2*i-1,2*j-1]<-0
          D2M[2*i,2*j]<-0
          D2M[2*i-1,2*j]<-CMax[i,i]
          D2M[2*j,2*i-1]<-D2M[2*i-1,2*j]
        }
      }
      else {
        D2M[2*i-1,2*j-1]<-CMin[i,j]
        D2M[2*i,2*j]<-CMax[i,j]
        D2M[2*i-1,2*j]<-(CMax[i,j]+CMin[i,j])/2
        D2M[2*i,2*j-1]<-(CMax[i,j]+CMin[i,j])/2
        D2M[2*j-1,2*i-1]<-D2M[2*i-1,2*j-1]
        D2M[2*j,2*i]<-D2M[2*i,2*j]
        D2M[2*j,2*i-1]<-D2M[2*i-1,2*j]
        D2M[2*j-1,2*i]<-D2M[2*i,2*j-1]        
      }
    }
  }
  return(D2M)
}
