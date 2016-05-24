BootFactorScores <-
function(PartialFS,niter = 1000){
# Bootstrap Confidence interval for DISTATIS
# computed on the partial factor scores
# PartialFS is nI,nF,nK array of the partial factor scores
# with nI # of objects, nF # of factors, nK # of observations
# (obtained from DISTATIS program)
# niter: how many iterations? default =1000
print(c('Bootstrap On Factor Scores. Iterations #: ',niter),quote=FALSE)


# Initialize the Bootstrap Table
nI <- dim(PartialFS)[1]
nK <- dim(PartialFS)[3]
nF <- dim(PartialFS)[2]
BootF <- array(0,dim = c(nI,nF,niter))
# Iterate Bootstrap
for (n in 1:niter){
	BootF[,,n] <- apply(PartialFS[,,sample(nK,nK,TRUE)],c(1,2),mean)
    }
rownames(BootF) <- rownames(PartialFS)
colnames(BootF) <- colnames(PartialFS)    
return(BootF)    
}
