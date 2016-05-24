nonzero.glmnet.cr <-
function(fit,s) {
	a0<-fit$a0[s]
	beta<-fit$beta[,s]
	beta<-beta[beta!=0]
	list(a0=a0,beta=beta)
}

