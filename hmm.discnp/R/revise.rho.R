revise.rho <- function(y,gamma,rnms) {
	y <- factor(unlist(y),levels=rnms)
	there <- !is.na(y)
	t1 <- apply(gamma[,there],1,
		function(x,index){tapply(x,index,sum)},y[there])
	t1[is.na(t1)] <- 0
	Rho <- t(t(t1)/apply(t1,2,sum))
        rownames(Rho) <- rnms
        Rho
}
