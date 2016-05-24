du.bcr <-
function(eta, Ymat, k, Cum.Ymat, x, link) {
	if (link=="logit") {
		z1 <- 1 / (1+exp(eta)) * Ymat[,2:k] + 
			exp(eta) / (1+exp(eta)) * (Cum.Ymat-1)
		#  swap in some C code?
		u <- t(x) %*% z1
        update.value <- min(rowSums(-u))
        update.j <- which.min(rowSums(-u))  
    } else if (link=="probit") {
        z1<- dnorm(eta)/pnorm(eta)*Ymat[,2:k] + dnorm(eta)/(1-pnorm(eta)+1e-16)*(Cum.Ymat-1)
		u<- -rowSums(t(x)%*%z1)
        update.value<-min(u)	
        update.j<-which.min(u)
	} else if (link=="cloglog") {
		d.delta<-matrix(0,ncol=k-1,nrow=dim(x)[1])
		for (i in 1:(k-1)) {
			d.delta[,i]<-exp(-exp(eta[,i])+eta[,i])
		}
		z1<- (1/G(eta, link))*Ymat[,2:k]*d.delta + (1-Cum.Ymat)*(1/(1-G(eta, link)))*-d.delta	
		u <- t(x) %*% z1
        update.value <- min(rowSums(-u, na.rm=TRUE))
        update.j <- which.min(rowSums(-u, na.rm=TRUE))  
	}
    list(update.j=update.j, update.value=update.value)  
}
