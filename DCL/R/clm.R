## CHAIN LADDER METHOD: parameter estimation with development factors and predictions
clm<-function(triangle,n.cal=NA,Fj=NA)
{
  Fj<-as.vector(Fj)
  n.cal<-as.vector(n.cal)
  triangle<-as.matrix(triangle)
  m<-nrow(triangle)
  # the arguments are 'triangle' the data with dimension 'm' 
  Z<-t(apply(triangle,1,cumsum)) # apply returns always the traspose
  # Ri are the sums by rows and Ci are the sums by columns
  Ri<-rowSums(triangle,na.rm=T)
  # Now if n.cal=NA and Fj is not provided or Fj=NA or length(Fj)<m-1 
  # we calculate forward factors as normally in chain ladder
  if (length(Fj)<m-1) {
    if (is.na(n.cal) ) {
      Cj<-colSums(triangle,na.rm=T)
      SRi<-cumsum(Ri)
      SCj<-cumsum(Cj[m:1])
      Si<-SRi-SCj
      # Development factors: Fj  (length=m-1 and we start in j=1,...,m-1
      Fj<-(Si[(m-1):1]+Cj[-1])/Si[(m-1):1]
      Fj[is.na(Fj)]<-1
      } else { 
          # if Fj is not a valid vector and n.cal is provided we calculate forward factors
          # Using only the most recent n.cal number of years  
          Top<-numeric(m-1)
          Bottom<-numeric(m-1)
          for (j in 1:(m-1))
          {
            i0<-max((m-j-n.cal+1),1)
            Top[j]<- sum(Z[i0:(m-j),j+1],na.rm=TRUE)
            Bottom[j]<-sum(Z[i0:(m-j),j],na.rm=TRUE)
          }
          # Developmen factors: Fj  (length=m-1 and we start in j=1,...,m-1
          Fj<-Top/Bottom
          Fj[is.na(Fj)]<-1
      } 
    } #else the provided development factors will be used in chain ladder
  
  beta<-numeric(m)  
  beta[1]<-1/prod(Fj,na.rm=TRUE)
  den<-cumprod(Fj[(m-1):1])[(m-1):1]
  beta[2:m]<-(Fj-1)/den
  ## The predicted lower and upper triangle
  Ypredict<-triangle
  # lower triangle
  Ypredict[2,m]<- Ri[2]*(Fj[(m-1)]-1)
  for (i in 3:m)
  {
    Ypredict[i,m-i+2]<- Ri[i]*(Fj[(m-i+1)]-1)
    for (j in (m-i+3):m)
    {
      Ypredict[i,j]<- Ri[i]*(Fj[(j-1)]-1) * prod(Fj[(m+1-i):(j-2)])
    }
  }
  #uppper triangle
  for (i in 1:m)
  { 
    # j=1
    if (i==m)  Ypredict[i,1]<- Ri[i] else Ypredict[i,1]<- Ri[i]* prod(1/Fj[1:(m-i)])
    if (i<m) {
      for (j in 2:(m-i+1)) Ypredict[i,j]<- Ri[i]*(Fj[(j-1)]-1) *prod(1/Fj[(j-1):(m-i)])  
    }
  }
  triangle.hat<-as.matrix(Ypredict)
  
  # now the alphas from the cummulative triangle.hat
  Z.hat<-t(apply(triangle.hat,1,cumsum)) ## apply returns always the traspose
  # Ri are the sums by rows and Ci are the sums by columns
  alpha<-Z.hat[,m]  
  par.clm<-list(triangle.hat=triangle.hat,alpha=alpha,beta=beta,Fj=Fj)
  
  return(par.clm)
}
