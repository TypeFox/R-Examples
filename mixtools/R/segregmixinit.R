segregmix.init <- function (y, x, lambda = NULL, beta = NULL, s = NULL, k = 2, seg.Z, 
psi, psi.locs = NULL) 
{
	n <- length(y)
	p <- ncol(x)
	psi.counts <- apply(psi>0,1,sum)
if(is.null(lambda)|is.null(beta)|is.null(s)|is.null(psi.locs)){
	if (is.null(psi.locs)) {
		psi.locs = vector("list",k)
		psi.locs = lapply(1:k, function(i) if(psi.counts[i]>0) vector("list",psi.counts[i]) else NULL)
		for(i in 1:k){
			if(!is.null(psi.locs[[i]])){
				temp.locs <- which(psi[i,]>0)
				temp.labs=NULL
				for(j in 1:length(temp.locs)){
					psi.locs[[i]][[j]]=sort(runif(psi[i,temp.locs[j]],as.numeric(quantile(x[,temp.locs[j]],.05)),as.numeric(quantile(x[,temp.locs[j]],.95))))
					temp.labs=c(temp.labs,colnames(x)[temp.locs[j]])
				}
			names(psi.locs[[i]])=temp.labs
			}
		}
	} else k = length(psi.locs)
	xnam <- colnames(x)
	fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+")))
	TEMP.lm <- lm(fmla,data=x)
	EM.out <- regmixEM(TEMP.lm$res,TEMP.lm$fit,k=k,epsilon=1e-2)
	posts = apply(EM.out$post,1,which.max)
	if (is.null(lambda)) {
		lambda = EM.out$lambda
		if(length(unique(posts))!=k) posts=rep(1:k,n)[1:n]
	} else k = length(lambda)
	A = round(lambda * n)
	while (min(A) <= 4) {
		lambda = runif(k)
		lambda = lambda/sum(lambda)
		A = round(lambda * n)
	}
	w = cbind(y, x)
	w.bin = list()
	for (j in 1:k) {
		w.bin[[j]] <- w[posts==j, ]
	}
	all.w.bin=vector("list",gamma(k+1))
	all.inds=perm(k,k)
	all.X.aug=all.w.bin
	all.lm.out=all.w.bin
	avg.res=NULL
	for(j in 1:length(all.w.bin)){
		all.w.bin[[j]]=w.bin[all.inds[j,]]
		X.aug <- lapply(1:k, function(i) cbind(1,aug.x(w.bin[[all.inds[j,i]]][,-1],unlist(psi.locs[[i]]),psi[i,],delta=NULL)))
		sapply(X.aug,dim)
		lm.out <- lapply(1:k, function(i) lm(w.bin[[all.inds[j,i]]][, 1] ~ X.aug[[i]][,-1]))
		all.X.aug[[j]]=X.aug
		all.lm.out[[j]]=lm.out
		avg.res=c(avg.res,mean(as.vector(unlist(lapply(1:k,function(t) lm.out[[t]]$res)))^2))
	}
	IND=which.min(avg.res)
	w.bin=all.w.bin[[IND]]
	X.aug=all.X.aug[[IND]]
	lm.out=all.lm.out[[IND]]
	s.hyp = lapply(lm.out, anova)
	s.hyp = as.vector(sqrt(sapply(1:k, function(i) tail(s.hyp[[i]]$Mean,1))))
	s.hyp[(s.hyp <= 0) | (is.na(s.hyp) == 1)] = 1
	if (is.null(s)) {
		s = 1/rexp(k, rate = s.hyp)
	} else k = length(s)
	if (is.null(beta)) {
	x.x <- lapply(1:k,function(i) try(solve(t(X.aug[[i]]) %*% X.aug[[i]]),silent=TRUE))
	test <- sum(sapply(1:k, function(i) class(x.x[[i]])[1]=="try-error"))
	if(test>0) stop("Lapack Routine Error")
	beta.hyp = lapply(lm.out,coef) # matrix(sapply(lm.out, coef), ncol = k)
	beta = vector("list",k)
	for (j in 1:k) {
		beta[[j]] = rnorm(length(beta.hyp[[j]]),mean=as.vector(beta.hyp[[j]]),
			sd = (s.hyp[j] * sqrt(diag(x.x[[j]]))))
	}
	} else k = length(beta)
} else{
		for(i in 1:k){
			if(!is.null(psi.locs[[i]])){
				temp.locs <- which(psi[i,]>0)
				temp.labs=NULL
				for(j in 1:length(temp.locs)){
					temp.labs=c(temp.labs,colnames(x)[temp.locs[j]])
				}
			names(psi.locs[[i]])=temp.labs
			}
		}
}
	list(lambda = lambda, beta = beta, s = s, k = k, psi.locs = psi.locs)
}

