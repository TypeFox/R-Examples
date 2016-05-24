distatis <-
function(LeCube2Distance, Norm = 'MFA', Distance = TRUE, 
                     RV = TRUE, nfact2keep = 3,compact = FALSE){
# DISTATIS
# Implement the standard DISTATIS program 
# described in Abdi H. et al, 2005, 2007, 2009, & 2012.
# (References to be completed)
# The compact option is used for the bootstrap

# Private functions
# DblCenterDist Create the Centering Matrix
DblCenterDist <- function(Y){ # Double Center a distance matrix
	nI = nrow(Y)
    CentMat = diag(nI) - (1/nI)*matrix(1,nI,nI)
	S = -(1/2)*(CentMat%*%Y%*%CentMat)
	return(S)
} # end of DblCenterDist
#
# *************************************************************************************************
Dist2CP <- function(D3){ # Transform a Cube of Distance into a cube of CP
CP3 <- (array(apply(D3,3,DblCenterDist),dim=c(dim(D3)[1],dim(D3)[2],dim(D3)[3])))
dimnames(CP3) <- dimnames(D3)
return(CP3)
} # end of Dist2CP
#
# *************************************************************************************************
# NormMFA4CP Normalize a CP matrix product to first eigenvalue of one
MFAnormCP <- function(Y){# MFA Normalize a psd matrix (i.e., first eigenvalue=1)
	ev = eigen(Y,symmetric=T,only.values=T)[1]
	e1 = ev$values[1]
	Ynormed = Y/e1
	return(Ynormed)
} # End of MFAnormed
#
# *************************************************************************************************
# CP2MFAnormedCP: MFA normalize a cube of CP
  CP2MFAnormedCP<- function(CP3){ # Transform a cube of CP into an MFA normed cube of CP
  CP3normed <- array(apply(CP3,3,MFAnormCP),dim=dim(CP3))
  dimnames(CP3normed) <- dimnames(CP3)
  return(CP3normed)
}
# *************************************************************************************************
# Compute the RV coefficient matrix
#
GetCmat <- function(CubeCP,RV=TRUE){# Compute a C matrix (or RV) from a cube of CP
	# get the dimensions
	nI = dim(CubeCP)[1]
	nJ = dim(CubeCP)[3]
	# reshape the 3D array as an (nI^2)*(nJ) matrix
	CP2 =  array(CubeCP,dim=c(nI*nI,nJ))
	C = t(CP2)%*%CP2 # Scalar Product
	if (RV){ # RV is TRUE we want the RV coefficient instead of the Scalar Product
		   laNorm = sqrt(apply(CP2^2,2,sum))
		    C = C / ( t(t(laNorm))%*%laNorm)
        } # end if
    rownames(C) <- dimnames(CubeCP)[[3]] -> colnames(C)      
    return(C)
}
#
# *************************************************************************************************
# Compute the compromise
ComputeSplus <- function(CubeCP,alpha){# Compute the comprise matrix for STATIS/DISTATIS
	nJ = dim(CubeCP)[3];
	nI = dim(CubeCP)[1]
	Splus = matrix(0,nI,nI)
	# Former Horrible Loop. Changed in final version 
	# for (i in 1:nJ){Splus = Splus + alpha[i]*CubeCP[,,i]}
	# A better way
	Splus <- apply( apply(CubeCP,c(1,2),'*',t(alpha)),c(2,3),sum)
	return(Splus)
    } # end ComputeSplus
# End of Private Routines
# *************************************************************************************************    
  # Transform the cube of distances into a cube of Cross-Products
  # if Distance is already a covariance or correlation reverse the sign
  if  (Distance != TRUE) {LeCube2Distance = -LeCube2Distance}
  # double center
  CP3 = Dist2CP(LeCube2Distance)
  # perform MFA normalization
  if (Norm == 'MFA') {CP3 <- CP2MFAnormedCP(CP3)}
  # Maybe more options here in the future)
  # Compute the C matrix as an RV matrix
  #
  C = GetCmat(CP3,RV=RV) # (to match matlab implementation of distatis) 
                         # get the RV coefficient mat
                         # instead of the scalar product 
  
  # eigen of C
  eigC <- eigen(C,symmetric=TRUE) # Eigen of C
  
  if (compact==FALSE){
     eigC$vectors[,1] <- abs(eigC$vectors[,1])
     rownames(eigC$vectors) <- rownames(C) 
     nfC = ncol(eigC$vectors) # number of eigenvectors of C
     Nom2Dim = paste('dim',1:nfC)
     colnames(eigC$vectors) <- Nom2Dim
     # make sure that first eigenvector is positive
     eigC$tau <- round(100*(eigC$values / sum(eigC$values)))
     names(eigC$tau)    <- Nom2Dim
     names(eigC$values) <- Nom2Dim
     # factors scores for RV mat
     eigC$G = t(apply(eigC$vectors,1,'*',t(t(sqrt(abs(eigC$values))))))
     rownames(eigC$G) <- rownames(C)
     colnames(eigC$G) <- Nom2Dim
  }
  # alpha weights
  alpha = eigC$vectors[,1] / sum(eigC$vectors[,1])
  # compute compromise
  Splus = ComputeSplus(CP3,alpha)
  if (compact == FALSE){
     # Eigen decomposition of Splus
     eigenSplus = eigen(Splus,symmetric=TRUE)
     # Percent Of Inertia
     eigenSplus$tau = round(100*eigenSplus$values / sum(eigenSplus$values))
     # singular values
     eigenSplus$SingularValues = sqrt(abs(eigenSplus$values))
     # Get the factor scores (ugly way, but works)
     F <- t(apply(eigenSplus$vectors,1,'*',t(t(eigenSplus$SingularValues))))
     rownames(F) <- rownames(Splus)
     # Keep only the interesting factors
     Nom2Factors = paste('Factor',1:nfact2keep)
     F <- F[,1:nfact2keep]
     colnames(F) <- Nom2Factors
     # Projection matrix
     ProjMat <- t(apply(eigenSplus$vectors,1,'*',1/t(t(eigenSplus$SingularValues))))
     Proj = ProjMat[,1:nfact2keep]
     colnames(Proj) <- Nom2Factors
     rownames(Proj) <- rownames(Splus)
    # Get the partial projections
    PartialF = array(apply(CP3,3,'%*%',Proj),
       dim=c(dim(CP3)[1],nfact2keep,dim(CP3)[3])    )
    rownames(PartialF) <- rownames(Splus)
    colnames(PartialF) <- Nom2Factors
    dimnames(PartialF)[[3]] <- rownames(C)     
       
    # pack up the information to send back
    # Will try to Keep a structure similar to  Cherise Chin Fatt MExPosition
    # in the meantime go for fast and match matlab
     	 res.Cmat <- list(C = C, eigVector = eigC$vector,eigValues = eigC$values,
                   tau = eigC$tau, G = eigC$G, alpha=alpha, compact=compact) 
         class(res.Cmat) <- c("Cmat","list")
    	 res.Splus <- list(SCP = CP3,F = F, PartialF = PartialF, ProjectionMatrix = Proj,Splus=Splus, compact=compact)
         class(res.Splus) <- c("Splus","list")    	 
		res.distatis <- list(res4Cmat= res.Cmat,res4Splus =res.Splus, compact=compact)
		class(res.distatis) <- c("DistatisR","list")
     } # End of if compact == FALSE
     else {# What do you do when it TRUE send back only the compact information 
    	 res.Cmat <- list(alpha=alpha, compact=compact) 
         class(res.Cmat) <- c("Cmat","list")    	 
	     res.Splus <- list(Splus=Splus, compact=compact)
         class(res.Splus) <- c("Splus","list")   	     
		 res.distatis <- list(res4Cmat= res.Cmat,res4Splus =res.Splus, compact=compact)	     
		class(res.distatis) <- c("DistatisR","list")		 
     }
     
     return(res.distatis)     
}
