#########################################################
#### AUTHOR:     Arnost Komarek                      ####
####             (2005)                              ####
####                                                 ####
#### FILE:       sampleCovMat.R                      ####
####                                                 ####
#### FUNCTIONS:  sampleCovMat                        ####
#########################################################

### ======================================
### sampleCovMat
### ======================================
sampleCovMat <- function(sample)
{
  dimnames <- colnames(sample)
  sample <- as.matrix(sample)

  n.sample <- nrow(sample)

  ybar <- apply(sample, 2, mean)              ## mean
  t.y.min.ybar <- t(sample) - ybar            ## (sample - mean)'
  covmat <- t.y.min.ybar %*% t(t.y.min.ybar)  ## sum (y-ybar)*(y-ybar)'

  if (n.sample > 1)  covmat <- (1/(n.sample - 1)) * covmat
  ## else return only sum of squares, which is zero

  rownames(covmat) <- dimnames
  colnames(covmat) <- dimnames  
  
  return(covmat)  
}  
