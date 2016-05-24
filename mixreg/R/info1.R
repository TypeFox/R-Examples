info1 <- function(x,y,theta,gma) {
#
# Note: The matrix x includes the column of 1's corresponding
# to the intercept term if an intercept is being fitted.
#
tmp <- list()
K   <- length(theta)
cK  <- sum(gma[,K])/theta[[K]]$lambda**2
po  <- function(a,b) {
	array(apply(b,2,"*",a), dim=c(dim(a),dim(b)[2]))
	}
tmp <- list()
for(k in 1:K) {
	bk <- theta[[k]]$beta
	vk <- theta[[k]]$sigsq
	lk <- theta[[k]]$lambda
	rk <- drop(y - x%*%bk)
# beta-beta:
	tmp1 <- -apply(gma[,k]*po(x,x),c(2,3),sum)/vk
# sigsq-sigsq:
	tmp2 <- -as.matrix(0.5*sum(gma[,k])/vk**2)
# lambda-lambda:
	tmp3 <- if(k < K) -as.matrix((sum(gma[,k])/lk**2 + cK)) else NULL
	tmp[[k]] <- dir.sum(tmp1,tmp2,tmp3)
}
#
# Note all ``cross terms'' are 0; beta-lambda and sigsq-lambda obviously
# so; sigsq-beta because the expression is essentially the sum of
# the weights times x times the residuals which is (linear algebra) 0.
#
rslt <- dir.sum(tmp)

ind <- (ncol(x)+2)*(1:(K-1))
n   <- length(ind)
m   <- cbind(rep(ind,n),rep(ind,rep(n,n)))
m   <- m[m[,1]!=m[,2],]
rslt[m] <- -cK

-rslt
}
