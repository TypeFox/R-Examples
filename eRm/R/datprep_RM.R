`datprep_RM` <-
function(X,W,sum0)                       #prepares data matrix for Rasch model
{ 
  X01 <- X                                        #X is already X(0,1)
  
  mt_vek <- rep(1,dim(X01)[2])                    #number of categories for each item
  K <- length(mt_vek)
  
  #automatized generation of the design matrix W
  if (length(W)==1) {
    W1 <- diag(1,(K-1))                           #build up design matrix
    if (sum0) {
      w1 <- rep(-1,(K-1))                         #sum0 restriction
    } else {
      w1 <- rep(0,(K-1))                          #first item parameter set to 0
    }
    W <- rbind(w1,W1)                             #RM design matrix  
    colnames(W) <- NULL
    rownames(W) <- NULL  
  }                                                     
  list(X=X,X01=X01,mt_vek=mt_vek,W=W)
#Output: X01      ... 0/1 response matrix of dimension N*rtot
#        mt_vek   ... 1-vector of length K 
#        W        ... design matrix of dimension K*K 
}

