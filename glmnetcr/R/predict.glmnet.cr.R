predict.glmnet.cr <-
function(object,newx=NULL, ...) {
	y <- object$y
	if (is.null(newx)) newx<-object$x
	method <- object$method
	j <- length(unique(y))  ### k in predict.glmpath.cr
	if (class(newx)=="numeric") newx<-matrix(newx,ncol=dim(object$x)[2])
	n <- dim(newx)[1]
	p <- dim(newx)[2]
	y.mat<-matrix(0,nrow=n,ncol=j)
	for(i in 1:n) y.mat[i,y[i]]<-1
	beta.est<-object$beta
	k<-apply(beta.est,2,function(x) sum(x!=0))   ### not in predict.glmpath.cr
	glmnet.BIC<-glmnet.AIC<-numeric()
	pi <-array(NA, dim=c(n,j,dim(beta.est)[2]))
	p.class<-matrix(NA,nrow=n,ncol=dim(beta.est)[2])
	for (i in 1:dim(beta.est)[2]) {
		beta<-beta.est[,i]
		logit <- matrix(rep(0, n * (j - 1)), ncol = j - 1)
		for (h in 1:(j - 1)) {
			cp <- paste("cp", h, sep = "")
			logit[, h] <- object$a0[i] + beta[names(beta) == cp] + beta[1:p] %*% t(as.matrix(newx))
		}
		delta <- matrix(rep(0, n * (j - 1)), ncol = j - 1)
		for (h in 1:(j - 1)) {
			delta[, h] <- exp(logit[, h])/(1 + exp(logit[, h]))
		}
		minus.delta <- 1 - delta
		if (c("backward", "forward")[charmatch(method, c("backward", 
														 "forward"))] == "backward") {
			for (h in j:2) {
				if (h == j) {
					pi[,h,i] <- delta[, j - 1]
				}
				else if (class(minus.delta[, h:(j - 1)]) == "numeric") {
					pi[,h,i] <- delta[, h - 1] * minus.delta[, h]
				}
				else if (dim(minus.delta[, h:(j - 1)])[2] >= 2) {
					pi[,h,i] <- delta[, h - 1] * apply(minus.delta[, 
													   h:(j - 1)], 1, prod)
				}
			}
			if (n==1) pi[, 1,i] <- 1 - sum(pi[, 2:j,i]) else pi[, 1,i] <- 1 - apply(pi[, 2:j,i], 1, sum)
		}
		if (c("backward", "forward")[charmatch(method, c("backward", 
														 "forward"))] == "forward") {
			for (h in 1:(j - 1)) {
				if (h == 1) {
					pi[,h,i] <- delta[, h]
				}
				else if (h == 2) {
					pi[,h,i] <- delta[, h] * minus.delta[, h - 1]
				}
				else if (h > 2 && h < j) {
					pi[,h,i] <- delta[, h] * apply(minus.delta[, 1:(h - 
																	1)], 1, prod)
				}
			}
			if (n==1) pi[, j,i] <- 1 - sum(pi[, 1:(j - 1),i]) else pi[, j,i] <- 1 - apply(pi[, 1:(j - 1),i], 1, sum)
		}
		LL<-0
		if (c("backward", "forward")[charmatch(method, c("backward", 
														 "forward"))] == "backward") {
			for (h in 1:(j-1)) {
				if (class(y.mat[,1:h])=="matrix") ylth<-apply(y.mat[,1:h],1,sum) else ylth<-y.mat[,1]
				LL<- LL + log(delta[,h])*y.mat[,h+1]+log(1-delta[,h])*ylth
			}
		}
		if (c("backward", "forward")[charmatch(method, c("backward", 
														 "forward"))] == "forward") {
			for (h in 1:(j-1)) {
				if (class(y.mat[,h:j])=="matrix") ygeh<-apply(y.mat[,h:j],1,sum) else ygeh<-y.mat[,j]
				LL<- LL + log(delta[,h])*y.mat[,h]+log(1-delta[,h])*ygeh
			}
		}
		LL<-sum(LL)
		glmnet.BIC[i]<- -2*LL+k[i]*log(n)
		glmnet.AIC[i]<- -2*LL+2*k[i]
		if (n==1) p.class[,i]<-which.max(pi[,,i]) else p.class[,i]<-apply(pi[,,i],1,which.max)
	}
	class<-matrix(levels(object$y)[p.class],ncol=ncol(p.class))
	names(glmnet.BIC)<-names(glmnet.AIC)<-dimnames(p.class)[[2]]<-dimnames(pi)[[3]]<-names(object$a0)
	dimnames(pi)[[2]]<-levels(object$y)
	list(BIC=glmnet.BIC,AIC=glmnet.AIC,class=class,probs=pi)		
}