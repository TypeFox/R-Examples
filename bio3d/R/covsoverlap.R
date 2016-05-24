covsoverlap <- function(...)
  UseMethod("covsoverlap")

covsoverlap.enma <- function(enma, ncore=NULL, ...) {
  if(!inherits(enma, "enma"))
    stop("provide a 'enma' object as obtain from function 'nma.pdbs()'")
  if(any(is.na(enma$fluctuations)))
    stop("provide 'enma' object calculated with argument 'rm.gaps=TRUE'")
  
  ncore <- setup.ncore(ncore, bigmem = FALSE)

  if(ncore>1)
    mylapply <- mclapply
  else
    mylapply <- lapply

  cat("Calculating pairwise covariance overlap coefs")

  m <- dim(enma$U.subspace)[3]
  mat <- matrix(NA, m, m)
  ##inds <- pairwise(m)
  inds <- rbind(pairwise(m),
                matrix(rep(1:m,each=2), ncol=2, byrow=T))
  
  mylist <- mylapply(1:nrow(inds), function(row) {
    i <- inds[row,1]; j <- inds[row,2];
    a <- list(U=enma$U.subspace[,,i], L=enma$L[i, ])
    b <- list(U=enma$U.subspace[,,j], L=enma$L[j, ])
    val <- covsoverlap.nma(a, b, ...)
    out <- list(val=val, i=i, j=j)
    cat(".")
    return(out)
  })
  
  for ( i in 1:length(mylist)) {
    tmp <- mylist[[i]]
    mat[tmp$i, tmp$j] <- tmp$val
  }

  mat[ inds[,c(2,1)] ] = mat[ inds ]
  ##diag(mat) <- rep(1, n)

  rownames(mat) <- basename(rownames(enma$fluctuations))
  colnames(mat) <- basename(rownames(enma$fluctuations))
  
  cat("\n")
  return(round(mat, 6))
}


covsoverlap.nma <- function(a, b, subset=NULL, ...) {
  if(any(missing(a), missing(b)))
    stop("provide eigenvectors and eigenvalues")

  dims.a <- dim(a$U)
  dims.b <- dim(b$U)
  
  if(dims.a[1]!=dims.b[1])
    stop("dimension mismatch")

  if(!is.null(subset)) {
    if(subset>ncol(a$U))
      subset <- ncol(a$U)
    
    a$U <- a$U[,1:subset]
    b$U <- b$U[,1:subset]
    a$L <- a$L[1:subset]
    b$L <- b$L[1:subset]
  }

  sumb <- 0
  for( k in 1:ncol(a$U) ) {
    tmp <- sqrt(a$L[k] * b$L)
    overlap <- c((t(a$U[,k]) %*% b$U)**2)
    sumb <- sumb + sum( tmp * overlap )
  }
  
  return(1 - ( sum(a$L + b$L) - 2 *sumb ) / sum(a$L + b$L))
}

