"mle" <-
function(v, r, ni)
{
# v is the var/cov matrix of paired differences
# r initially should be the grubbs estimate of the 
#   squared precision (vector 1 to N)
# r when returned is the MLE
# N is the no. of measuring devices
# ni is the no. of measured items

  N <- length(r)

  for(i in 1:N)
  {
    b0 <- 0
    b1 <- 0
    b2 <- 0
  
    for(j in 1:N)
    {
      if(i != j)
      {
        b0 <- b0 + 1/r[j]
        b1 <- b1 + v[i,j]/r[j]
        jplus1 <- j + 1
        if(jplus1 <= N)
        {
          for(k in jplus1:N)
          {
            if(k != i) b2 <- b2 + v[j,k]/(r[j]*r[k])
          }
        }
      }
    }
    r[i] <- (((ni-1)*(b0*b1-b2))/(ni*b0^2))-1/b0
#   cat("i = ",i, " r = ",r[i], " b0 = ", b0, " b1 = ",b1, " b2 = ",b2,"\n",sep="")
  }
  r
}
