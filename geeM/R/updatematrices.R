### Get the AR-1 inverse matrix
getAlphaInvAR <- function(alpha.new, a1,a2,a3,a4, row.vec, col.vec){
  corr.vec <- c(alpha.new*a1/(1-alpha.new^2) , ( (1+alpha.new^2)*a2 + a3)/(1-alpha.new^2) + a4, alpha.new*a1/(1-alpha.new^2))	
  return(as(sparseMatrix(i= row.vec, j=col.vec, x=corr.vec), "symmetricMatrix"))
}


### Get the necessary structures for the inverse matrix of AR-1
### returns the row and column indices of the AR-1 inverse
### a1, a2, a3, and a4 are used to compute the entries in the matrix
buildAlphaInvAR <- function(len){
  
  nn <- sum(len)
  K <- length(len)
  a1 <- a2 <- a3 <- a4 <- vector("numeric", nn)
  index <- c(cumsum(len) - len, nn)
  
  for (i in 1:K)  {
    if(len[i] > 1)  {
      a1[(index[i]+1) : index[i+1]] <- c(rep(-1,times = len[i]-1),0)
      a2[(index[i]+1) : index[i+1]] <- c(0,rep(1,times=len[i]-2),0)
      a3[(index[i]+1) : index[i+1]] <- c(1,rep(0,times=len[i]-2),1)
      a4[(index[i]+1) : index[i+1]] <- c(rep(0,times=len[i]))
    }
    else if (len[i] == 1)  {
      a1[(index[i]+1) : index[i+1]] <- 0
      a2[(index[i]+1) : index[i+1]] <- 0
      a3[(index[i]+1) : index[i+1]] <- 0
      a4[(index[i]+1) : index[i+1]] <- 1
    }
  }
  
  a1 <- a1[1:(nn-1)]
  subdiag.col <- 1:(nn-1)
  subdiag.row <- 2:nn
  
  row.vec <- c(subdiag.row, (1:nn), subdiag.col)
  col.vec <- c(subdiag.col, (1:nn), subdiag.row)
  return(list(row.vec = row.vec, col.vec= col.vec, a1 = a1, a2=a2, a3=a3, a4=a4))
}

### Returns the full inverse matrix of the correlation for EXCHANGEABLE structure
getAlphaInvEX <- function(alpha.new, diag.vec, BlockDiag){
  return(as(BlockDiag %*% Diagonal(x = (-alpha.new/((1-alpha.new)*(1+(diag.vec-1)*alpha.new)))) + Diagonal( x = ((1+(diag.vec-2)*alpha.new)/((1-alpha.new)*(1+(diag.vec-1)*alpha.new)) + alpha.new/((1-alpha.new)*(1+(diag.vec-1)*alpha.new)))), "symmetricMatrix"))
}

### Get the inverse of M-DEPENDENT correlation matrix.
getAlphaInvMDEP <- function(alpha.new, len, row.vec, col.vec){
  K <- length(len)
  N <- sum(len)
  m <- length(alpha.new)
  
  # First get all of the unique block sizes.
  mat.sizes <- sort(unique(len))
  corr.vec <- vector("numeric", sum(len^2))
  mat.inverses <- list()
  index <- c(0, (cumsum(len^2) -len^2)[2:K], sum(len^2)) 
  
  for(i in 1:length(mat.sizes)){
    # Now create and invert each matrix
    if(mat.sizes[i] == 1){
      mat.inverses[[i]] <- 1
    }else{		
      mtmp <- min(m, mat.sizes[i]-1)
      a1 = list()
      a1[[1]] <- rep(1, mat.sizes[i])
      for(j in 1:mtmp){
        a1[[j+1]] <- rep(alpha.new[j], mat.sizes[i]-j)
      }
      
      tmp <- bandSparse(mat.sizes[i], k=c(0:mtmp),diagonals=a1, symmetric=T )
      mat.inverses[[i]] <- as.vector(solve(tmp))
    }
  }
  # Put all inverted matrices in a vector in the right order
  corr.vec <- unlist(mat.inverses[len - min(len) + 1])
  return(as(sparseMatrix(i=row.vec, j=col.vec, x=corr.vec), "symmetricMatrix"))
}

### Get a vector of elements of the inverse correlation matrix for UNSTRUCTURED
### Inversion strategy follows the same basic idea as M-DEPENDENT
getAlphaInvUnstruc <- function(alpha.new, len, row.vec, col.vec){
  K <- length(len)
  unstr.row <- NULL
  unstr.col <- NULL
  ml <- max(len)
  sl2 <- sum(len^2)
  for(i in 2:ml){
    unstr.row <- c(unstr.row, 1:(i-1))
    unstr.col <- c(unstr.col, rep(i, each=i-1))
  }
  unstr.row <- c(unstr.row, 1:ml)
  unstr.col <- c(unstr.col, 1:ml)
  xvec <- c(alpha.new, rep(1, ml))
  # Get the biggest matrix implied by the cluster sizes
  biggestMat <- forceSymmetric(sparseMatrix(i=unstr.row, j=unstr.col, x=xvec))
  
  mat.sizes <- sort(unique(len))
  corr.vec <- vector("numeric", sl2)
  mat.inverses <- list()
  index <- vector("numeric", K+1)
  index[1] <- 0
  index[2:K] <-  (cumsum(len^2) -len^2)[2:K]
  index[K+1] <- sl2
  
  for(i in 1:length(mat.sizes)){
    tmp <- biggestMat[1:mat.sizes[i], 1:mat.sizes[i]]
    mat.inverses[[i]] <- as.vector(solve(tmp))		
  }
  
  corr.vec <- unlist(mat.inverses[len - min(len) + 1])
  return(as(sparseMatrix(i=row.vec, j=col.vec, x=corr.vec), "symmetricMatrix"))
  
}



### Invert the FIXED correlation structure.  Again,
### uses same basic technique as M-DEPENDENT
getAlphaInvFixed <- function(mat, len){
  K <- length(len)
  mat.sizes <- sort(unique(len))
  mat.inverses <- list()
  sl2 <- sum(len^2)	
  corr.vec <- vector("numeric", sl2)
  index <- vector("numeric", K+1)
  index[1] <- 0
  index[2:K] <-  (cumsum(len^2) -len^2)[2:K]
  index[K+1] <- sl2
  for(i in 1:length(mat.sizes)){
    tmp <- mat[1:mat.sizes[i], 1:mat.sizes[i]]
    mat.inverses[[i]] <- as.vector(solve(tmp))
  }
  
  corr.vec <- unlist(mat.inverses[len - min(len) + 1])
  
  return(as(getBlockDiag(len, corr.vec)$BDiag, "symmetricMatrix"))	
}

