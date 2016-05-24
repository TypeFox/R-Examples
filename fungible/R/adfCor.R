##########################################################
# This function compute the asymptotic distribution-free #
# covariance matrix of correlations.                     #
#                                                        #
# Arguments:                                             #
# X - matrix of predictor scores                         #
# y - vector of criterion scores                         #
#                                                        #
# Output                                                 #
# adfCovMat - Asymptotic distribution-free estimate of   #
#             the covariance matrix                      #
##########################################################

adfCor <- function(X,y=NULL) {
	
######################## Internal Functions ########################	
	
# row.remove: Removes rows from a symmetric transition matrix 
# to create a correlation transition matrix  
# (see Browne & Shapiro, (1986); Nel, 1985).

  row.remove <- function(p) {
    p1 <- p2 <- p
    rows <- rep(1,p)
    for(i in 2:p) {
      rows[i] <- rows[i] + p1
      p1 <- p1 + (p2-1)
      p2 <- p2 - 1
    }
    rows
  }
  
# Dp: Duplicator Matrix	 
  Dp <- function(p) {
    M <- matrix(nrow = p, ncol = p)
    M[ lower.tri(M, diag = T) ] <- seq( p*(p + 1)/2 )
    M[ upper.tri(M, diag = F) ] <- t(M)[ upper.tri(M, diag = F) ]
    D <- outer(c(M), unique(c(M)),
                FUN = function(x, y) as.numeric(x == y) )
    D
  }  

###################### End Internal Functions ######################	
####################################################################	
  Xy <- if(is.null(y)) X else cbind(X,y)
  dev <- scale(Xy, scale=FALSE)
  nvar <- ncol(dev)
  N <- nrow(dev)

# order of the covariance matrix of covariances
  ue <- nvar^2

# container for indices
  s <- vector(length=ue, mode="character")

  z <- 0
  for(i in 1:nvar){
    for(j in 1:nvar){
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
  ids <- matrix(0,nrow=ue^2,4)
  for(i in 1:4) ids[,i] <- as.numeric(sapply(V,substr,i,i))
  
  covs <- matrix(0,(ue^2),1)

# compute the covariances using Steiger and Hakstian (1982) Eqn 3.4
  for(i in 1:(ue^2)) {

    w_ii <- cov(dev[,ids[i,1]],dev[,ids[i,1]])
    w_jj <- cov(dev[,ids[i,2]],dev[,ids[i,2]])
    w_kk <- cov(dev[,ids[i,3]],dev[,ids[i,3]])
    w_ll <- cov(dev[,ids[i,4]],dev[,ids[i,4]])

    w_ij <- cov(dev[,ids[i,1]],dev[,ids[i,2]])
    w_kl <- cov(dev[,ids[i,3]],dev[,ids[i,4]])
    w_ik <- cov(dev[,ids[i,1]],dev[,ids[i,3]])
    w_jl <- cov(dev[,ids[i,2]],dev[,ids[i,4]])
        
    w_ijkl <- t(dev[,ids[i,1]]*dev[,ids[i,2]])%*%(dev[,ids[i,3]]*dev[,ids[i,4]])/(N-1)
    w_iikk <- t(dev[,ids[i,1]]*dev[,ids[i,1]])%*%(dev[,ids[i,3]]*dev[,ids[i,3]])/(N-1)
    w_jjkk <- t(dev[,ids[i,2]]*dev[,ids[i,2]])%*%(dev[,ids[i,3]]*dev[,ids[i,3]])/(N-1)
    w_iill <- t(dev[,ids[i,1]]*dev[,ids[i,1]])%*%(dev[,ids[i,4]]*dev[,ids[i,4]])/(N-1)
    w_jjll <- t(dev[,ids[i,2]]*dev[,ids[i,2]])%*%(dev[,ids[i,4]]*dev[,ids[i,4]])/(N-1)
    w_iikl <- t(dev[,ids[i,1]]*dev[,ids[i,1]])%*%(dev[,ids[i,3]]*dev[,ids[i,4]])/(N-1)    
    w_jjkl <- t(dev[,ids[i,2]]*dev[,ids[i,2]])%*%(dev[,ids[i,3]]*dev[,ids[i,4]])/(N-1)    
    w_ijkk <- t(dev[,ids[i,1]]*dev[,ids[i,2]])%*%(dev[,ids[i,3]]*dev[,ids[i,3]])/(N-1)    
    w_ijll <- t(dev[,ids[i,1]]*dev[,ids[i,2]])%*%(dev[,ids[i,4]]*dev[,ids[i,4]])/(N-1)    
            
    r_ij <- w_ij/sqrt(w_ii*w_jj)
    r_kl <- w_kl/sqrt(w_kk*w_ll)
    r_ik <- w_ik/sqrt(w_ii*w_kk)
    r_jl <- w_jl/sqrt(w_jj*w_ll) 
    
    r_ijkl <- w_ijkl/sqrt(w_ii*w_jj*w_kk*w_ll)
    r_iikk <- w_iikk/sqrt(w_ii*w_ii*w_kk*w_kk)
    r_jjkk <- w_jjkk/sqrt(w_jj*w_jj*w_kk*w_kk)           
    r_iill <- w_iill/sqrt(w_ii*w_ii*w_ll*w_ll)
    r_jjll <- w_jjll/sqrt(w_jj*w_jj*w_ll*w_ll)
    r_iikl <- w_iikl/sqrt(w_ii*w_ii*w_kk*w_ll)
    r_jjkl <- w_jjkl/sqrt(w_jj*w_jj*w_kk*w_ll)                                                        
    r_ijkk <- w_ijkk/sqrt(w_ii*w_jj*w_kk*w_kk)           
    r_ijll <- w_ijll/sqrt(w_ii*w_jj*w_ll*w_ll)               
        
    covs[i] <- (r_ijkl + .25*r_ij*r_kl*(r_iikk + r_jjkk + r_iill + r_jjll) -
               .5*r_ij*(r_iikl + r_jjkl) - .5*r_kl*(r_ijkk + r_ijll))
  }

  adfFull <- matrix(covs,ue,ue,byrow=T)  

# Create symmetric transition matrix (see Nel, 1985)
  Kp <- solve(t(Dp(nvar)) %*% Dp(nvar)) %*% t(Dp(nvar))  
  
# Create correlation transition matrix (see Browne & Shapiro, 1986).
  Kpc <- as.matrix(Kp[-row.remove(nvar),] )
  if(ncol(Kpc) == 1) Kpc <- t(Kpc)  
  
  adfCovMat <- Kpc %*% adfFull %*% t(Kpc)/N
  
  ## add row and column labels
  rc<-expand.grid(1:nvar,1:nvar)
  rc<-rc[!(rc[,1]==rc[,2]),]
  rc<-unique(apply(t(apply(rc,1,sort)),1,paste, collapse=""))
  dimnames(adfCovMat)<-list(rc,rc)
  adfCovMat

}