pcr <-
function(x, nTreat, M, m){

  if(!is.data.frame(x))              stop('x must be a data.frame')
  if(nTreat < 1 | nTreat >= nrow(x)) stop('nTreat is number to be treated: 
                                          need 0 < nTreat < nrow(x)')
  if(M < 1)                          stop('M must be a positive integer 
                                          (suggested M >= 10,000)')
  if(m < 1 | m >= M)                 stop('need 0 < m < M 
                                          (suggested m/M <= 0.10)')

  tIndexMat = genPerms(n=nrow(x), n1=nTreat, nPerms=M)
  vars = sapply(1:ncol(tIndexMat), function(i){
    getVar(covs=x, tIndex=tIndexMat[,i])
  })
  
  mCut = quantile(vars, m/M)
  mKeep = which(vars <= mCut)

  return(list(treatments=tIndexMat, variance=vars, cutoff=mCut, best=mKeep))
}
