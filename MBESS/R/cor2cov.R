`cor2cov` <-
function(cor.mat, sd, discrepancy=1e-5)
{
if(dim(cor.mat)[1]!=dim(cor.mat)[2]) stop("'cor.mat' should be a square matrix")

n<-sqrt(length(cor.mat))
if(n!=length(sd)) stop("The length of 'sd' should be the same as the number of rows of 'cor.mat'")

if(length(sd[sd>0])!= n) stop("The elements in 'sd' shuold all be non-negative")

if(isSymmetric(cor.mat)) IS.symmetric <- TRUE 
	else IS.symmetric <- FALSE
p <- dim(cor.mat)[1]
q <- p*(p-1)/2
if (isTRUE(all.equal(cor.mat[lower.tri(cor.mat)], rep(0,q))) || isTRUE(all.equal(cor.mat[upper.tri(cor.mat)], rep(0,q)))) 
	IS.triangular <- TRUE
	else IS.triangular <- FALSE
if(!IS.symmetric & !IS.triangular)	stop("The object 'cor.mat' should be either a symmetric or a triangular matrix")

cov.mat <- diag(sd)  %*% cor.mat  %*% diag(sd)
colnames(cov.mat)<- rownames(cov.mat)<- colnames(cor.mat)
return(cov.mat)
}