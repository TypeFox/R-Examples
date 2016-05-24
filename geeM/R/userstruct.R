
### Get the structure of the USERDEFINED correlation matrix implied by the
### corr.mat argument to geem.
getUserStructure <- function(corr.mat){
  ml <- dim(corr.mat)[1]
  
  row.vec <- NULL
  col.vec <- NULL
  for(i in 2:ml){
    row.vec <- c(row.vec, 1:(i-1))
    col.vec <- c(col.vec, rep(i, each=i-1))
  }
  
  struct.vec <- corr.mat[cbind(row.vec, col.vec)]
  
  corr.list <- vector("list", max(struct.vec))
  for(i in 1:max(struct.vec)){
    corr.list[[i]] <- which(struct.vec == i)
  }
  return(list(corr.list = corr.list, row.vec = row.vec, col.vec = col.vec, struct.vec = struct.vec))
}


### Get the inverse correlation matrix for USERDEFINED.
getAlphaInvUser <- function(alpha.new, len, struct.vec, user.row, user.col, row.vec, col.vec){
  K <- length(len)
  ml <- max(len)
  sl2 <- sum(len^2)
  
  # Indices for the correlation matrix for the subject
  # with the most observations.
  user.row <- c(user.row, 1:ml)
  user.col <- c(user.col, 1:ml)
  # The entries of the biggest matrix
  xvec <- rep.int(0, length(struct.vec))
  for(i in 1:length(alpha.new)){
    xvec[which(struct.vec == i)] <- alpha.new[i]
  }
  
  xvec <- c(xvec, rep(1, ml))	
  
  biggestMat <- forceSymmetric(sparseMatrix(i=user.row, j=user.col, x=xvec))
  
  mat.sizes <- sort(unique(len))
  corr.vec <- vector("numeric", sl2)
  mat.inverses <- list()
  
  for(i in 1:length(mat.sizes)){
    tmp <- biggestMat[1:mat.sizes[i], 1:mat.sizes[i]]
    mat.inverses[[i]] <- as.vector(solve(tmp))		
  }
  
  
  corr.vec <- unlist(mat.inverses[len - min(len) + 1])
  return(as(sparseMatrix(i=row.vec, j=col.vec, x=corr.vec), "symmetricMatrix"))
}

### Check some conditions on the USERDEFINED correlation structure supplied.
checkUserMat <- function(corr.mat, len){
  if(is.null(corr.mat)){
    stop("corr.mat must be specified if using user defined correlation structure")
  }
  if(dim(corr.mat)[1] < max(len)){
    stop("corr.mat needs to be at least as long as the maximum cluster size.")
  }
  test.vec <- as.vector(corr.mat)
  if(any(abs(test.vec-round(test.vec)) > .Machine$double.eps )){
    stop("entries in corr.mat must be integers.")
  }
  max.val <- max(test.vec)
  min.val <- min(test.vec)
  if(!all(sort(unique(test.vec)) == min.val:max.val)){
    stop("entries in corr.mat must be consecutive integers starting at 1.")
  }
  return(corr.mat[1:max(len), 1:max(len)])	
}
