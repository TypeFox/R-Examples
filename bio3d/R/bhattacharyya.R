bhattacharyya <- function(...)
  UseMethod("bhattacharyya")

bhattacharyya.nma <- function(...)
  bhattacharyya.matrix(...)

bhattacharyya.pca <- function(...)
  bhattacharyya.matrix(...)

bhattacharyya.enma <- function(enma, covs=NULL, ncore=NULL, ...) {
  if(!inherits(enma, "enma"))
    stop("provide a 'enma' object as obtain from function 'nma.pdbs()'")
  if(any(is.na(enma$fluctuations)))
    stop("provide 'enma' object calculated with argument 'rm.gaps=TRUE'")
  
  if(is.null(covs)) {
    cat("Calculating covariance matrices")
    covs <- cov.enma(enma, ncore=ncore)
  }
  cat("Calculating pairwise bhattacharyya coefs")
  sim.mat <- bhattacharyya.array(covs, ncore=ncore)
  
  rownames(sim.mat) <- basename(rownames(enma$fluctuations))
  colnames(sim.mat) <- basename(rownames(enma$fluctuations))
  return(sim.mat)
}

bhattacharyya.matrix <- function(a, b, q=90, n=NULL, ...) {

  if(!is.matrix(a) & is.matrix(b))
    stop("provide covariance matrices")
  
  a <- ((1 / .tr(a)) * a)*1000
  b <- ((1 / .tr(b)) * b)*1000
  ei <- eigen( (a + b)/2  )

  if(is.null(n)) {
    percent <- (ei$values/sum(ei$values))*100
    cumv <- cumsum(percent)
    n <- which(cumv>q)[1]
  }

  ca <- det((t(ei$vectors[,1:n]) %*% a) %*% ei$vectors[,1:n])
  cb <- det((t(ei$vectors[,1:n]) %*% b) %*% ei$vectors[,1:n])

  d <- prod(ei$values[1:n])
  
  ndb <- (1/(2*n)) * log( d / sqrt(ca*cb) )
  bc <- exp( -ndb )
  return(bc)
}

bhattacharyya.array <- function(covs, ncore=NULL, ...) {
  ncore <- setup.ncore(ncore, bigmem = FALSE)

  if(ncore>1)
    mylapply <- mclapply
  else
    mylapply <- lapply

  dims <- dim(covs)
  m <- dims[3]
  
  mat <- matrix(NA, m, m)
  ##inds <- pairwise(m)
  inds <- rbind(pairwise(m),
                matrix(rep(1:m,each=2), ncol=2, byrow=T))
                
  mylist <- mylapply(1:nrow(inds), function(row) {
    i <- inds[row,1]; j <- inds[row,2];
    val <- bhattacharyya.matrix(covs[,,i], covs[,,j], ...)
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
  
  cat("\n")
  return(round(mat, 6))
}
