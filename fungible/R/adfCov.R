##########################################################
# This function compute the asymptotic distribution-free #
# covariance matrix of covariances.                      #
#                                                        #
# Arguments:                                             #
# X - matrix of predictor scores                         #
# y - vector of criterion scores                         #
#                                                        #
# Output                                                 #
# adfCovMat - Asymptotic distribution-free estimate of   #
#             the covariance matrix                      #
##########################################################

adfCov <- function(X, y=NULL) {
  Xy <- if(is.null(y)) X else cbind(X,y)
  dev <- scale(Xy,scale=FALSE)
  nvar <- ncol(dev)
  N <- nrow(dev)

# order of the covariance matrix of covariances
  ue <- nvar*(nvar + 1)/2

# container for indices
  s <- vector(length=ue, mode="character")

  z <- 0
  for(i in 1:nvar){
    for(j in i:nvar){
      z<-z+1
      s[z]<-paste(i,j,sep="")
    }
  }

# computes all possible combinations of the 
# indices in s
  v <- expand.grid(s, s)

# paste the index pairs togehter
  V <- paste(v[,1], v[,2], sep="")

# separate the indices into their own columns
  id.mat <- matrix(0,nrow=ue^2,4)
  for(i in 1:4) id.mat[,i] <- as.numeric(sapply(V,substr,i,i))

# create a matrix with the sequence of numbers 1:ue^2 by row
# we will use this to find the positions of the indices in id.mat
  M <- matrix(1:ue^2, ue, ue, byrow=TRUE)

# grabs the rows of the index pairs of interest
  r <- M[lower.tri(M, diag=TRUE)]

  ids <- id.mat[r,]

  adfCovMat <- adfCovBias <- matrix(0,ue,ue)
  covs <- covs.bias <- matrix(0,nrow(ids),1)

# compute the covariances using Browne (1984) Eqn 3.8

  for(i in 1:nrow(ids)) {
   
    w_ij <- crossprod(dev[,ids[i,1]],dev[,ids[i,2]])/N
    w_ik <- crossprod(dev[,ids[i,1]],dev[,ids[i,3]])/N
    w_il <- crossprod(dev[,ids[i,1]],dev[,ids[i,4]])/N
    w_jk <- crossprod(dev[,ids[i,2]],dev[,ids[i,3]])/N
    w_jl <- crossprod(dev[,ids[i,2]],dev[,ids[i,4]])/N
    w_kl <- crossprod(dev[,ids[i,3]],dev[,ids[i,4]])/N
     
  
    w_ijkl <- (t(dev[,ids[i,1]]*dev[,ids[i,2]])%*%
                 (dev[,ids[i,3]]*dev[,ids[i,4]])/N)         
    
    covs[i] <- (N*(N-1)*(1/((N-2)*(N-3)))*(w_ijkl - w_ij*w_kl) -
                  N*(1/((N-2)*(N-3)))*(w_ik*w_jl + w_il*w_jk - (2/(N-1))*w_ij*w_kl))
  }	

# Create Covariance Matrix

  adfCovMat[lower.tri(adfCovMat,diag=T)] <- covs	
  vars <- diag(adfCovMat)
  adfCovMat <- adfCovMat + t(adfCovMat) - diag(vars)

  ## add row and column labels
  rc<-expand.grid(1:nvar,1:nvar)
  rc<-unique(apply(t(apply(rc,1,sort)),1,paste, collapse=""))
  rc
  dimnames(adfCovMat)<-list(rc,rc)
  adfCovMat
  
}
