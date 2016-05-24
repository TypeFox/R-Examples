testrightHAND<-function(U)
{
  ## require(RFOC)
  K = matrix(c(1,2,3,2,3,1,3,1,2), ncol=3, byrow=TRUE)
  tes = rep(NA, 3)
for(i in 1:3)
  {
    

    x12  = RSEIS::xprod(U[  ,K[i,1] ] , U[,K[i,2] ])
      x3 = U[,K[i,3]]
      tes[i] = all(sign(x3)==sign(x12))
      
  }
return(tes)
  
}
