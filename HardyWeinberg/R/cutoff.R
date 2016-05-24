cutoff <-
function(n,alpha=0.05,cc=0.5,verbose=FALSE) {
  chiquant <- qchisq(1-alpha,1)
  zposnum <- -1.5*cc^2 + cc*sqrt((9*(cc^2))/4 + 4*n*chiquant)
  zposden <- 2*n*chiquant
  zpos <- zposnum/zposden
  if(verbose) cat("zpos = ",zpos,"\n")
  root <- 0.5-0.5*sqrt(1-4*zpos)
#  cc2 <- cc*cc
#  k1 <- cc2
#  k2 <- -1.5*cc2
#  k3 <- -n*chiquant
#  coef <- c(k1,k2,k3)
#  out <- polyroot(coef)
#  print(out)
  return(root)
}

