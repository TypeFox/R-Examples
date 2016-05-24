`datprep_LPCM` <-
function(X,W,mpoints,Groups,sum0)
{
  #TFrow <- (rowSums(X)==0)                       #el. persons with 0 rawscore
  #X <- X[!TFrow,]

  ngroups <- max(Groups)
  N <- dim(X)[1]                                  #number of persons
  K <- dim(X)[2]/mpoints                          #number of items
  mt_vek <- apply(X,2,max,na.rm=TRUE)[1:K]                   #number of categories - 1 for each item
  mt_vek_0 <- mt_vek+1                            #number of categories for each item
  
  X01_0 <- matrix(rep(0,(N*sum(mt_vek_0)*mpoints)),nrow=N) #empty 0/1 matrix  
  K1 <- dim(X)[2]
  cummt0 <- c(0,cumsum(rep(mt_vek_0,mpoints))[1:(K1-1)])+1     #index vector for 0th category
  indmatp <- apply(X,1,function(xi) {xi+cummt0})  #preparing index matrix for 1 responses
  imp1 <- as.vector(indmatp)
  imp2 <- rep(1:N,rep(K1,N))
  indmat <- cbind(imp2,imp1)                      #final index matrix for 1 responses
  X01_0[indmat] <- 1                              #0/1 matrix with 0th category
  
  
  d1 <- 1:N
  d2 <- 1:K1
  coor <- expand.grid(d2,d1)[,c(2:1)]               #X coordinates
  resvec <- as.vector(t(X))                         #X as vector (rowwise)
  NAind <- as.matrix(coor[is.na(resvec),])          #index matrix for NA's in X
  mt_vek.t <- rep(mt_vek,mpoints)
    
  if (length(NAind) > 0) {
    NAindlist <- apply(NAind,1,function(x){
                    #x <- unlist(x)
                    co <- seq(cummt0[x[2]],cummt0[x[2]]+mt_vek.t[x[2]])
                    NAind01 <- cbind(rep(x[1],length(co)),co)
                    rownames(NAind01) <- NULL
                    data.frame(NAind01,row.names=NULL)                                               #list with NA indices
                    })
    indmatNA <- matrix(unlist(lapply(NAindlist, function(x) {t(as.matrix(x))})),ncol=2,byrow=TRUE)   #matrix with NA indices 
    X01_0[indmatNA] <- NA
  }
  
  X01 <- X01_0[,-cummt0]
    
  #automatized generation of the design matrix W
  if (length(W)==1) {
    W11diag <- diag(1,(sum(mt_vek)-1))                   #build up design matrix
    
    if (sum0) {
      w110 <- rep(-1,(sum(mt_vek)-1))                 #sum0 restriction
    } else {
      w110 <- rep(0,(sum(mt_vek)-1))                  #first item category parameter set to 0
    }
    
    W11 <- rbind(w110,W11diag)                               #PCM design matrix 
    ZW <- dim(W11)[1]
    
    W1 <- NULL
    for (i in 1:(mpoints*ngroups)) W1 <- rbind(W1,W11)  #first part with virtual items
    
    if (mpoints > 1) {                            #more than 1 measurement points
      if (ngroups > 1) {                          #more than 1 group/more mpoints
        t_mp1 <- rep(1:mpoints,rep(ZW*ngroups,mpoints))
        t_mp <- factor(t_mp1)
        g_ng1 <- rep(rep(1:ngroups,rep(ZW,ngroups)),mpoints)
        g_ng <- factor(g_ng1)
        W2 <- model.matrix(~t_mp+g_ng)[,-1]               #main effects g and mp
        W2[1:(ZW*ngroups),] <- 0                          #remove main effects for the first test occasion 
      } else {                                    #1 group/more mpoints
        t_mp <- gl(mpoints,ZW)             #factor for measurement points
        W2 <- model.matrix(~t_mp)[,-1] }
    } else if (ngroups > 1) {                     #1 mpoint/more groups
        g_ng <- gl(ngroups,ZW)
        W2 <- model.matrix(~g_ng)[,-1] 
        warning("Group contrasts without repeated measures can not be estimated!")
    } else if (ngroups == 1) W2 <- NULL           #1 mpoint/1 group
        
  catvek <- sequence(mt_vek)
  W2_cat <- W2*catvek                             #imposing item categories
  W <- cbind(W1,W2_cat)                           #design matrix completed
  colnames(W) <- NULL
  rownames(W) <- NULL 
  }
  
   
  list(X=X,X01=X01,mt_vek=mt_vek,W=W)
}

