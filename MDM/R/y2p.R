y2p <- function(y, mean = FALSE) {
	dn <- dimnames(y)
	y  <- as.matrix(y)/apply(y,1,sum)
	if (mean) y <- matrix(rep(apply(y,2,mean),each=nrow(y)),nrow=nrow(y),ncol=ncol(y),dimnames=dimnames(y))
	dimnames(y) <- dn
	y
}
