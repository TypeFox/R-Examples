fillout <-
function(z) {
	if(is.vector(z)) z <- matrix(z,ncol=1)
	qq <- nrow(z)
	l <- ncol(z)
	if(l<qq) {
		qz <- diag(qq)-z%*%solve(t(z)%*%z,t(z))
		z <- cbind(z,eigen(qz)$vector[,1:(qq-l)])
	}
	if(l>qq) {
		z <- t(fillout(t(z)))
	}
	z
}
