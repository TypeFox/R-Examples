semimetric.pca <-
function(DATA1, DATA2, q)
{
	if(is.vector(DATA1)) DATA1 <- as.matrix(t(DATA1))
	if(is.vector(DATA2)) DATA2 <- as.matrix(t(DATA2))
	testfordim <- sum(dim(DATA1)==dim(DATA2))==2
	twodatasets <- T
	if(testfordim) twodatasets <- sum(DATA1==DATA2)!=prod(dim(DATA1))
	qmax <- ncol(DATA1)
	if(q > qmax) stop(paste("give a integer q smaller than ", qmax))
	n <- nrow(DATA1)
	COVARIANCE <- t(DATA1) %*% DATA1/n
	EIGENVECTORS <- eigen(COVARIANCE, symmetric = TRUE)$vectors[, 1:q]
	COMPONENT1 <- DATA1 %*% EIGENVECTORS
	if(twodatasets) {
		COMPONENT2 <- DATA2 %*% EIGENVECTORS
	}
	else {
		COMPONENT2 <- COMPONENT1
	}
	SEMIMETRIC <- 0
	for(qq in 1:q)
		SEMIMETRIC <- SEMIMETRIC + outer(COMPONENT1[, qq], COMPONENT2[, 
			qq], "-")^2
	return(sqrt(SEMIMETRIC))
}
