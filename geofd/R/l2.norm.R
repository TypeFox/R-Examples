l2.norm <-
function(s, datafd, M)
{

  L2norm <- matrix(0,nrow=s,ncol=s)
  coef <- datafd$coef

  for (i in 1:(s-1)){
    coef.i <- coef[,i]
    for (j in (i+1):s){
      coef.j <- coef[,j]
      L2norm[i,j] <- t(coef.i-coef.j) %*% M %*% (coef.i-coef.j)
      L2norm[j,i] <- L2norm[i,j]
    }
  }

##################################################################
# Return:
##################################################################

  return(L2norm)

}
