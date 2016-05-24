BootFromCompromise <-
function(LeCube2Distance,niter =1000, Norm = 'MFA', 
     Distance = TRUE, RV = TRUE, nfact2keep = 3){
#  Bootstrap Confidence interval for DISTATIS
# computed on the partial factor scores
# PartialFS is nI,nF,nK array of the partial factor scores
# with nI # of objects, nF # of factors, nK # of observations
# (obtained from DISTATIS program)
# niter: how many iterations? default =1000
print(c('Starting Full Bootstrap. Iterations #: ',niter),quote=FALSE)

# First call distatis to get the fixed effect model
FixedDist <- distatis(LeCube2Distance, Norm=Norm, Distance=Distance,RV=RV,nfact2keep=nfact2keep)
# Projection Matrix    	
  ProjMat = FixedDist$res4Splus$ProjectionMatrix
 # Initialize the Bootstrap Table
nI <- dim(LeCube2Distance)[1]
nK <- dim(LeCube2Distance)[3]
nF <- min(c(dim(ProjMat)[2],nfact2keep))
FullBootF <- array(0,dim = c(nI,nF,niter))
rownames(FullBootF) <- rownames(LeCube2Distance)
colnames(FullBootF) <- paste('Factor',1:nF)
# Iterate Bootstrap
for (n in 1:niter){
	ResBoot_n  <- distatis(LeCube2Distance[,,sample(nK,nK,TRUE)], Norm=Norm,Distance=Distance,RV=RV,nfact2keep=nfact2keep,compact=TRUE)
	FullBootF[,,n] <- ResBoot_n$res4Splus$Splus %*% ProjMat
    }
return(FullBootF)  
}
