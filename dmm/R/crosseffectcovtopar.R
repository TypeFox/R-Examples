crosseffectcovtopar <-
function(covx, covt, varcomp, var1, var2){
# crosseffectcovtopar() - covs to parameters for cross-effect cases
#                       - same trait and cross trait cases
  n <- ncol(covx)  # n traits
  if(ncol(covt) != n) {
    stop("Covx and Covt must be same size in crosseffectcovtopar():\n")
  }
  if(ncol(varcomp) != n*n) {
    stop("Coviama and Varcomp must be same size in crosseffectcovtopar():\n")
  }
  fract <- rep(0,n)
  corre <- matrix(0,n,n)
  for(j in 1:n) {
    jj <- (j-1)*n + j
    for(i in 1:n) {
      ii <- (i-1)*n + i
      if(i == j) {
       fract[i] <- covx[i,i]/covt[i,i]
      }
      able <- varcomp[var1,jj] * varcomp[var2,ii]
      if(able > 0) {
        corre[i,j] <- covx[i,j]/sqrt(able)
      }
      else {
        corre[i,j] <- 0
      }
    }
  }
  return(list(corre=corre, fract=fract))
}
