`datprep_LLTM` <-
function(X,W,mpoints,Groups,sum0)
{
# Design matrix see Fischer & Molenaar, p. 159  
  
  #TFrow <- (rowSums(X)==0 | rowSums(X)==(dim(X)[2]))  #el. persons with 0/K rawscore
  #X <- X[!TFrow,]
  
  ngroups <- max(Groups)
  X01 <- X
  N <- dim(X)[1]                                  #number of persons
  K <- dim(X)[2]/mpoints                            #number of items
  mt_vek <- rep(1,K)
  
  #automatized generation of the design matrix W
  if (length(W)==1) {
    W11diag <- diag(1,(sum(mt_vek)-1))                #build up design matrix
    if (sum0) {
      w110 <- rep(-1,(sum(mt_vek)-1))                 #sum0 restriction
    } else {
      w110 <- rep(0,(sum(mt_vek)-1))                  #first item category parameter set to 0
    }
    W11 <- rbind(w110,W11diag)                        #RM design matrix 
    ZW <- dim(W11)[1]
    
    W1 <- NULL
    for (i in 1:(mpoints*ngroups)) W1 <- rbind(W1,W11)    #first part with virtual items
    
    if (mpoints > 1) {                                    #more than 1 measurement points
      if (ngroups > 1) {                                  #more than 1 group/more mpoints
        t_mp1 <- rep(1:mpoints,rep(ZW*ngroups,mpoints))
        t_mp <- factor(t_mp1)
        g_ng1 <- rep(rep(1:ngroups,rep(ZW,ngroups)),mpoints)
        g_ng <- factor(g_ng1)
        W2 <- model.matrix(~t_mp+g_ng)[,-1]               #main effects g and mp
        W2[1:(ZW*ngroups),] <- 0                          #remove main effects for the first test occasion 
      } else {                                            #1 group/more mpoints
        t_mp <- gl(mpoints,ZW)                            #factor for measurement points
        W2 <- model.matrix(~t_mp)[,-1] }
    } else if (ngroups > 1) {                             #1 mpoint/more groups
        g_ng <- gl(ngroups,ZW)
        W2 <- model.matrix(~g_ng)[,-1] 
        warning("Group contrasts without repeated measures can not be estimated!")
    } else if (ngroups == 1) W2 <- NULL                   #1 mpoint/1 group
        
  W <- cbind(W1,W2)
  colnames(W) <- NULL
  rownames(W) <- NULL 
  }
  
  list(X=X,X01=X01,mt_vek=mt_vek,W=W)
#Output: X01      ... 0/1 response matrix of dimension N*rtot
#        mt_vek   ... vector of length K with number of categories - 1 (for each item)
#        W        ... design matrix of dimension (K*T)*((K-1)*(T-1)+1)
}

