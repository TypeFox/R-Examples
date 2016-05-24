#########################################################
#### AUTHOR:     Arnost Komarek                      ####
####             (2005)                              ####
####                                                 ####
#### FILE:       vecr2matr.R                         ####
####                                                 ####
#### FUNCTIONS:  vecr2matr                           ####
#########################################################

### ======================================
### vecr2matr
### ======================================
### 26/01/2005
vecr2matr <- function(vec.r, KK)
{
  dim <- length(KK)
  if (dim < 1 || dim > 2) stop("Not implemented for dimension of the G-spline different from 1 or 2")

  if (dim == 1){
    init.r <- as.vector(vec.r) - 1 - KK[1]
    return(init.r)
  }
  if (dim == 2){
    r1 <- as.vector(vec.r - 1) %%  (2*KK[1] + 1) - KK[1]
    r2 <- as.vector(vec.r - 1) %/% (2*KK[1] + 1) - KK[2]
    result <- cbind(r1, r2)
    colnames(result) <- c("r1", "r2")
    return(result)
  }  
}  
