### subset the FPCA object with specified number of components k
### and return the subsetted fpcaObj

SubsetFPCA <- function(fpcaObj, k){
  fpcaObj$lambda <- fpcaObj$lambda[1:k]
  fpcaObj$phi <- fpcaObj$phi[,1:k, drop=FALSE]
  fpcaObj$xiEst <- fpcaObj$xiEst[,1:k, drop=FALSE]
  fpcaObj$FVE <- fpcaObj$cumFVE[k]
  fpcaObj$cumFVE <- fpcaObj$cumFVE[1:k]
  return(fpcaObj)
}
