IntermediateONN <-
function(plist, skew.vec, kurto.vec, ONNCorrMat) 
{
  cmat.corrected <- matrix(NA, nrow(ONNCorrMat), ncol(ONNCorrMat))
  
  coef.est <- Fleishman.coef.NN(skew.vec, kurto.vec)  
  # Compute the Fleishman coefficients for standard normal variable X
  coef.X <- Fleishman.coef.NN(0,0)
      
  Cor_forONN <- function(pvec, ONN.cor, coef) {
    X <- rnorm(1e+05, 0, 1)
    XORD <- ordinalize(pvec, X)
    c.hat <- cor(XORD, X) 
    rho.XY <- ONN.cor/c.hat
    
    # Given X is standard normal and Y is non-normal, need to compute the intermediate 
    # correlation between the two standard normal variables X and Z
    
    c1 <- -rho.XY
    c2 <- coef.X[2]*coef[2]+3*coef.X[2]*coef[4]+3*coef.X[4]*coef[2]+9*coef.X[4]*coef[4]
    c3 <- 2*coef.X[3]*coef[3]
    c4 <- 6*coef.X[4]*coef[4]
    roots <- polyroot(c(c1,c2,c3,c4))  
    for (k in 1:3) {
      if(abs(Im(roots[k]))<1e-6 & abs(Re(roots[k]))<=1)  {
        r <- Re(roots[k])
      }
    }      
    return(r)
  }
  
  # Correlation between the jth binary/ordinal variable and the ith non-normal variable
  
  for (j in 1:ncol(ONNCorrMat)) {
    for (i in 1:nrow(ONNCorrMat)) {
      cmat.corrected[i, j] = Cor_forONN(plist[[j]], ONNCorrMat[i,j], as.vector(coef.est[i,]))
    }
  }
  return(cmat.corrected)
}
