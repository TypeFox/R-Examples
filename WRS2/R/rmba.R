rmba<-function(x, csteps = 5,na.rm=TRUE,plotit=FALSE)
{
# computes the reweighted MBA estimator
# Code supplied by David Olive
#
#       x is assumed to be a matrix
#
#  plotit=FALSE is used to avoid problems when this function is called
#  by other function in WRS
#
x=as.matrix(x)
if(na.rm)x=elimna(x)
	p <- dim(x)[2]
	n <- dim(x)[1]	##get the DGK estimator
	covs <- var(x)
	mns <- apply(x, 2, mean)	## concentrate
	for(i in 1:csteps) {
		md2 <- mahalanobis(x, mns, covs)
		medd2 <- median(md2)
		mns <- apply(as.matrix(x[md2 <= medd2,  ]), 2, mean)
		covs <- var(x[md2 <= medd2,  ])
	}
	covb <- covs
	mnb <- mns	##get the square root of det(covb)
	critb <- prod(diag(chol(covb)))	##get the resistant estimator
	covv <- diag(p)
	med <- apply(x, 2, median)
	md2 <- mahalanobis(x, center = med, covv)
	medd2 <- median(md2)	## get the start
	mns <- apply(as.matrix(x[md2 <= medd2,  ]), 2, mean)
	covs <- var(x[md2 <= medd2,  ])	## concentrate
	for(i in 1:csteps) {
		md2 <- mahalanobis(x, mns, covs)
		medd2 <- median(md2)
#		mns <- apply(as.matrix(x[md2 <= medd2,  ]), 2, mean)
		mns <- apply(as.matrix(x[md2 <= medd2,  ]), 2, mean)
		covs <- var(x[md2 <= medd2,  ])
	}
	crit <- prod(diag(chol(covs)))
	if(crit < critb) {
		critb <- crit
		covb <- covs
		mnb <- mns
	}
##scale for better performance at MVN
	rd2 <- mahalanobis(x, mnb, covb)
	const <- median(rd2)/(qchisq(0.5, p))
	covb <- const * covb	
	##reweight the above MBA estimator (mnb,covb) for efficiency
	rd2 <- mahalanobis(x, mnb, covb)
	up <- qchisq(0.975, p)
	rmnb <- apply(as.matrix(x[rd2 <= up,  ]), 2, mean)
	rcovb <- var(x[rd2 <= up,  ])
	rd2 <- mahalanobis(x, rmnb, rcovb)
	const <- median(rd2)/(qchisq(0.5, p))
	rcovb <- const * rcovb	## reweight again
	rd2 <- mahalanobis(x, rmnb, rcovb)
	up <- qchisq(0.975, p)
	rmnb <- apply(as.matrix(x[rd2 <= up,  ]), 2, mean)
	rcovb <- var(x[rd2 <= up,  ])
	rd2 <- mahalanobis(x, rmnb, rcovb)
	const <- median(rd2)/(qchisq(0.5, p))
	rcovb <- const * rcovb
cor.b=NULL
temp=outer(sqrt(diag(rcovb)),sqrt(diag(rcovb)),'*')
if(min(diag(rcovb)>0))cor.b=rcovb/temp
	list(center = rmnb, cov = rcovb, cor=cor.b)
}