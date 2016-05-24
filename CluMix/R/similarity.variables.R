similarity.variables <-
function(data, associationFun=association, check.psd=TRUE, make.psd=TRUE){
# data: data.frame of original data 
# associationFun: function that calculates association measure for each pair of variables
# check.psd: check if resulting similarity matrix S is positive semi-definite?
# make.psd: if S is not p.s.d., shall it be transformed to be p.s.d.? (only done if also check.psd=TRUE)
  
  #n <- nrow(data)
  p <- ncol(data)
  
  S <- matrix(0, nrow=p, ncol=p)
  for(i in 1:p){
    for(j in 1:p){
      if(i > j){
        # distance = sqrt(1 - association)
        S[i,j] <- associationFun(data[,i], data[,j])
      }
    }
  }
  dimnames(S) <- list(names(data), names(data))
  
  # make it symmetric (since only lower triangle was calculated) 
  S <- S + t(S) 
  diag(S) <- 1
  
  # check if S is p.s.d.
  if(check.psd){
    psd <- all(eigen(S, only.values=TRUE)$values >= 0)
  
    # if S is not p.s.d., get "nearest p.s.d. matrix"
    if(!psd){
      if(!make.psd)
        warning("similarity matrix is not positive semidefinite")
      else{
        S <- Matrix::nearPD(S, keepDiag=TRUE, conv.norm.type="F")$mat
        #warning("similarity matrix was adjusted to be positive semidefinite")
      }
    }
  }
  
  return(S)
}
