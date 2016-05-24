intercor.NN <-
function(pmat, cmat){
  if (dim(pmat)[1] != dim(cmat)[1]) {
    stop("Correlation matrix dimension is not consistent with number of continous variables in pmat!\n")
  }
  
  if (dim(pmat)[2] != 4){
    stop("column of pmat must be 4\n")
  }
  
  cmat_star = diag(1,dim(cmat)[1])
  
  k = 1
  for (i in 2:dim(cmat)[1])
  {
    for (j in 1:(i-1)){
      z = rep(NA,3)
      z[1]=pmat[i,2]*pmat[j,2]+3*pmat[i,2]*pmat[j,4]+3*pmat[i,4]*pmat[j,2]+9*pmat[i,4]*pmat[j,4]
      z[2]=2*pmat[i,3]*pmat[j,3]
      z[3]=6*pmat[i,4]*pmat[j,4]
      
      c_star = polyroot(c(-cmat[i,j],z[1],z[2],z[3]))
      
      c_star = c_star[abs(Im(c_star))<1e-10 & abs(Re(c_star))<=1]
      cmat_star[i,j] = cmat_star[j,i] = Re(c_star)
      k = k+1
    }
  }
  return(round(cmat_star,3))
}
