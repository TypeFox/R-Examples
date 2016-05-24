cov.nma <- function(nma) {
  if(!inherits(nma, "nma"))
    stop("provide a 'nma' object as obtain from function 'nma.pdb()'")
  
  dims <- dim(nma$U)
  cov <- matrix(0, ncol=dims[1], nrow=dims[1])
  tmpU <- nma$U[, (nma$triv.modes+1):ncol(nma$U)]
  tmpL <- nma$L[(nma$triv.modes+1):ncol(nma$U)]
  
  for(j in 1:ncol(tmpU) ) {
    cov <- cov + ( (tmpU[,j] %*% t(tmpU[,j])) / tmpL[j])
  }
  return(cov)
}

cov.enma <- function(enma, ncore=NULL) {
  if(!inherits(enma, "enma"))
    stop("provide a 'enma' object as obtain from function 'nma.pdbs()'")
  if(any(is.na(enma$fluctuations)))
    stop("provide 'enma' object calculated with argument 'rm.gaps=TRUE'")
  
  ncore <- setup.ncore(ncore, bigmem = FALSE)
  
  if(ncore>1)
    mylapply <- mclapply
  else
    mylapply <- lapply
  
  if(!inherits(enma, "enma"))
    stop("provide 'enma' object as obtained from nma.pdbs")
  
  dims <- dim(enma$U.subspace)
  
  mycalc <- function(i, enma) {
    cov <- matrix(0, ncol=dims[1], nrow=dims[1])
    tmpU <- enma$U.subspace[,,i]
    tmpL <- enma$L[i,]

    for(j in 1:ncol(tmpU) ) {
      cov = cov + ( (tmpU[,j] %*% t(tmpU[,j])) / tmpL[j])
    }
    cat(".")
    return(cov)
  }

  covs.list <- mylapply(1:dims[3L], mycalc, enma)
  cat("\n")
  
  covs <- array(0, dim=c(dims[1], dims[1], dims[3]))
  
  for ( i in 1:dims[3L] )
    covs[,,i]=covs.list[[i]]

  return(covs)
}

.tr <- function(mat) {
  return(sum(diag(mat)))
}

