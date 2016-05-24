#' Simulate Z based on a probit model
#' 
#' Simulates a random latent matrix Z given its expectation, dyadic correlation
#' and a binary relational matrix Y
#' 
#' 
#' @usage rZ_bin_fc(Z, EZ, rho, Y)
#' @param Z a square matrix, the current value of Z
#' @param EZ expected value of Z
#' @param rho dyadic correlation
#' @param Y square binary relational matrix
#' @return a square matrix , the new value of Z
#' @author Peter Hoff
#' @export rZ_bin_fc
rZ_bin_fc <-
function(Z,EZ,rho,Y)
{
  # simulates Z under the contraints
  # (1)  Y[i,j]=1   => Z[i,j]>0    
  # (2)  Y[i,j]=0   => Z[i,j]<0



  sz<-sqrt(1-rho^2)
  ut<-upper.tri(EZ)
  lt<-lower.tri(EZ)
 
  Y[is.na(Y)]<- -1  
  for(y in c((-1):1))
  {
    lb<-c(-Inf,-Inf,0)[y+2] ; ub<-c(Inf,0,Inf)[y+2]  

    up<- ut & Y==y
    ez<- EZ[up] + rho*( t(Z)[up]  - t(EZ)[up] )
    Z[up]<-ez+sz*qnorm(runif(sum(up),pnorm((lb-ez)/sz),pnorm((ub-ez)/sz)))

    up<- lt & Y==y
    ez<- EZ[up] + rho*( t(Z)[up]  - t(EZ)[up] )
    Z[up]<-ez+sz*qnorm(runif(sum(up),pnorm((lb-ez)/sz),pnorm((ub-ez)/sz)))
  }

  diag(Z)<-rnorm(nrow(Z),diag(EZ),1)
  Z
}
