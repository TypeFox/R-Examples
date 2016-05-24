intercor.PP <-
function(lamvec, cmat) {
  if (length(lamvec) != dim(cmat)[1]) {
    stop("Correlation matrix dimension is not consistent with number of Poisson variables in lamvec!\n")
  }
  if (sum(lamvec <= 0) > 0) {
    stop("lambda should be positive \n")
  }
  norow = 1e+05
  u = runif(norow, 0, 1)
  
  corrected = diag(1,length(lamvec))
  for (i in 2:length(lamvec)){
    for (j in 1:(i-1)){

      maxcor = cor(qpois(u, lamvec[i]), qpois(u,lamvec[j]))
      mincor = cor(qpois(u, lamvec[i]), qpois(1 - u,lamvec[j]))
      
      a = -maxcor * mincor/(maxcor + mincor)
      b = log((maxcor + a)/a, exp(1))
      c = -a
      corrected[i,j] = corrected[j,i] = log((cmat[i,j] + a)/a, exp(1))/b
      corrected[i,j] = corrected[j,i] = ifelse((corrected[i,j] > 1 | corrected[i,j] < (-1)), 
                                               NA, corrected[i,j])
    }    
  }
  return(round(corrected,3))
}
