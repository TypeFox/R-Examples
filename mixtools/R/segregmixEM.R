segregmixEM=function (y, x, lambda = NULL, beta = NULL, sigma = NULL, k = 2, seg.Z, psi, psi.locs = NULL, delta = NULL,
	epsilon = 1e-08, maxit = 10000, verb = FALSE, max.restarts=15) 
{
	if (sum(x[,1]==1)==nrow(x)) x=x[,-1]
	x=data.frame(x) 
	col.names.x <- colnames(x) 
	xnam <- colnames(x) 
	fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+"))) 
if(!is.null(psi.locs)){
psi.counts=apply(psi,1,sum)
for(i in 1:k){
	if(psi.counts[i]>0){
		TEMP <- (is.list(psi.locs[[i]]))&(length(psi.locs[[i]])==sum(psi[i,]>0)) 
	} else{
		TEMP <- is.null(psi.locs[[i]])
		}
	if(TEMP==FALSE) stop(paste("You must specify a correct changepoint structure!", 
		"\n"))
}
}

	if (!is.null(delta)) {
		cat("Estimation performed assuming the changepoints are known.", "\n")
		if (is.null(psi.locs)) {
			stop(paste("You must specify the changepoints for this setting!", 
				"\n"))
		}
	} 
	if ((length(seg.Z) != k) | class(seg.Z) != "list") {
		stop(paste("You must specify a list of length k for the segmented relationships!", 
			"\n"))
	}
	if (!identical(all.equal(dim(psi),c(k,ncol(x))),TRUE)) {
		stop(paste("You must specify a matrix with the correct dimension for psi!", 
			"\n"))
	}
	if (((length(psi.locs) != k) | class(psi.locs) != "list") & !is.null(psi.locs)) {
		stop(paste("You must specify a list of length k for the number of changepoints per predictor in each component!", 
			"\n"))
	}
	tot.cp <- apply(psi,1,sum)
	tmp.ind=1
	tmp <- try(suppressWarnings(segregmix.init(y=y, x=x, lambda = lambda, beta = beta, s = sigma, k = k, seg.Z=seg.Z, 
		psi=psi, psi.locs = psi.locs)),silent=TRUE)
	if(class(tmp)=="try-error"){
		cat("Generating new initial values.", "\n")
		while(tmp.ind<=10){
			tmp <- try(suppressWarnings(segregmix.init(y=y, x=x, lambda = NULL, beta = NULL, s = NULL, k = k, seg.Z=seg.Z, 
				psi=psi, psi.locs = NULL)),silent=TRUE)
			tmp.ind <- tmp.ind+1
			if(tmp.ind==11)  stop(paste("Had to reinitialize algorithm too many times.  Reconsider specified model.", "\n"))
			if(class(tmp)!="try-error") tmp.ind=20
		}
	}
	x.old=x
	x = cbind(1, x)
	data.x=cbind(y,x) 
	lambda <- tmp$lambda
	beta <- tmp$beta
	s <- tmp$s
	k <- tmp$k
	psi.locs <- tmp$psi.locs
	sing <- 0
	perms=perm(k,k)
	perm.ind=nrow(perms)
	if (is.null(delta)) delta <- lapply(1:k,function(i) NULL)
	n <- length(y)
	diff <- 1
	iter <- 0
	X.aug <- lapply(1:k, function(i) cbind(1,aug.x(x[,-1],unlist(psi.locs[[i]]),psi[i,],delta=delta[[i]])))
	X.aug.old <- X.aug
	psi.locs.old <- psi.locs
	xbeta <- lapply(1:k, function(i) X.aug[[i]] %*% matrix(beta[[i]],ncol=1))
	res <- sapply(1:k, function(i) as.numeric((y - xbeta[[i]])^2))
	comp <- t((lambda/sqrt(2 * pi * s^2)) * t(exp(-t(t(res)/(2 * s^2)))))
	obsloglik <- sum(log(apply(comp, 1, sum)))
	ll <- obsloglik
	z = matrix(nrow = n, ncol = k)
	restarts <- 0
	while (diff > epsilon && iter < maxit) {
	null.beta=0
	lambda.old <- lambda
	beta.old <- beta
	s.old <- s
        for (i in 1:n) {
		for (j in 1:k) {
			z.denom = c()
			for (h in 1:k) {
				z.denom = c(z.denom, (lambda[h]/lambda[j]) * 
					(s[j]/s[h]) * exp(-0.5 * ((1/s[h]^2) * res[i, h] - (1/s[j]^2) * res[i, j])))
			}
			z[i, j] = 1/sum(z.denom)
		}
	}
	z = z/apply(z, 1, sum)
	z[,k] = 1-apply(as.matrix(z[,-k]),1,sum)
	z<-round(z,15)
	z.old <- z
	lambda.new <- apply(z, 2, mean)
	if (sum(lambda.new < 1e-08) > 0 || is.na(sum(lambda.new))) {
		sing <- 1
	} else {
		lm.out <- vector("list",k)
		psi.temp=psi.locs
		psi.ind=lapply(1:k,function(i) which(psi[i,]!=0))
		for(i in 1:k){
			if(is.null(seg.Z[[i]]) | (sum(1-sapply(delta,is.null))>0)){
				temp.seg <- lm(fmla,data=data.x,weights=z[,i])
			} else temp.seg <- try(suppressWarnings(segmented(lm(fmla,data=data.x,weights=z[,i]),seg.Z=seg.Z[[i]],psi=psi.temp[[i]])),silent=TRUE)
			if(class(temp.seg)[1]=="try-error"){
				seq = 1
				temp.names = names(psi.locs.old[[i]])
				while(seq < 20){
					psi.temp2 <- vector("list",length(psi.temp[[i]]))
					for(ii in 1:length(psi.temp[[i]])){
						x.range <- range(data.x[,which(names(data.x)==temp.names[ii])])
						psi.temp2[[ii]] <- psi.temp[[i]][[ii]]+sample(c(-1,1),length(psi.temp[[i]][[ii]]),replace=TRUE)*runif(length(psi.temp[[i]][[ii]]),0,diff(x.range)/10)
						if((any(psi.temp2[[ii]]<=x.range[1]))|(any(psi.temp2[[ii]]>=x.range[2]))) psi.temp2[[ii]]=psi.temp[[i]][[ii]]
						psi.temp2[[ii]]=sort(psi.temp2[[ii]])
					}
					names(psi.temp2)=temp.names
					temp.seg <- try(suppressWarnings(segmented(lm(fmla,data=data.x,weights=z[,i]),seg.Z=seg.Z[[i]],psi=psi.temp2[[i]],control=seg.control(it.max=1))),silent=TRUE)
					if(class(temp.seg)[1]=="try-error"){
						seq = seq+1
					} else seq=40
				} 
				if(seq!=40){
					temp.seg <- try(suppressWarnings(segmented(lm(fmla,data=data.x,weights=z[,i]),seg.Z=seg.Z[[i]],psi=psi.temp[[i]],control=seg.control(it.max=1))),silent=TRUE)
				}
			}
			lm.out[[i]]=temp.seg			
		}
		lambda <- lambda.new
		if(sum(sapply(lm.out,class)=="try-error")>0){
			newobsloglik=-Inf
		} else{
			if(sum(1-sapply(delta,is.null))>0){
			psi.new <- psi.locs.old
			} else {
				psi.new <- psi.locs
				for(i in 1:k){
					if(class(lm.out[[i]])[1]=="segmented"){
						temp.names=names(psi.locs[[i]])
						temp.cumsum=cumsum(sapply(psi.locs[[i]],length))
						TC.ind = length(temp.cumsum)
						seg.temp = lm.out[[i]]$psi[,2]
						psi.new[[i]] = lapply(1:length(psi.locs[[i]]), function(j) as.numeric(lm.out[[i]]$psi[,2]))
						psi.new[[i]] = vector("list",TC.ind)
						psi.new[[i]][[1]]=sort(seg.temp[1:temp.cumsum[1]])
						if(TC.ind>1) for(j in 2:TC.ind) psi.new[[i]][[j]] = sort(seg.temp[(temp.cumsum[j-1]+1):temp.cumsum[j]])
						names(psi.new[[i]])=temp.names
					}
				}
			}
			X.aug.new <- lapply(1:k, function(i) cbind(1,aug.x(x[,-1],unlist(psi.new[[i]]),psi[i,],delta[[i]])))
			lm.out2=lapply(1:perm.ind, function(j) lapply(1:k, function(i) lm(y~X.aug.new[[i]][,-1],weights=z[,perms[j,i]])))
			beta.new <- lapply(1:perm.ind, function(j) lapply(lm.out2[[j]],coef))
			null.perms <- sapply(1:perm.ind,function(i) all(!is.na(lapply(beta.new,unlist)[[i]])))
			null.beta=0
			if(sum(null.perms)>0){
				xbeta.new <- lapply(1:perm.ind, function(j) lapply(1:k, function(i) X.aug.new[[i]] %*% matrix(beta.new[[j]][[i]],ncol=1)))
				res <- lapply(1:perm.ind, function(j) sapply(1:k, function(i) (y - xbeta.new[[j]][[i]])^2))
				s.new <- lapply(1:perm.ind, function(j) sqrt(sapply(1:k, function(i) sum(z[,perms[j,i]] * (res[[j]][, i]))/sum(z[,perms[j,i]]))))
				comp <- lapply(1:perm.ind, function(j) lapply(1:k, function(i) lambda.new[i] * dnorm(y, xbeta.new[[j]][[i]], s.new[[j]][i])))
				comp <- lapply(1:perm.ind, function(j) sapply(comp[[j]], cbind))
				compsum <- lapply(1:perm.ind, function(j) apply(comp[[j]], 1, sum))
				newobsloglik <- sapply(1:perm.ind, function(j) sum(log(compsum[[j]])))
				newobsloglik[c(1-null.perms)]=-Inf
				IND <- which.max(newobsloglik)
				z = z[,perms[IND,]]
				lambda.new <- apply(z, 2, mean)
				lambda <- lambda.new
				beta <- beta.new[[IND]]
				xbeta <- xbeta.new[[IND]]
				res <- res[[IND]]
				X.aug <- X.aug.old
				s <- s.new[[IND]]
				psi.locs <- psi.new
				sing <- sum(s < 1e-08)
				newobsloglik <- newobsloglik[IND]
			} else{
				newobsloglik=Inf
				null.beta=1
			} 
		}
		if((newobsloglik<obsloglik & !is.na(newobsloglik))|null.beta==1){
		lm.out1=lapply(1:perm.ind, function(j) lapply(1:k, function(i) lm(y~X.aug.old[[i]][,-1],weights=z[,perms[j,i]])))
		beta.new <- lapply(1:perm.ind, function(j) lapply(lm.out1[[j]],coef))
		null.perms <- sapply(1:perm.ind,function(i) all(!is.na(lapply(beta.new,unlist)[[i]])))
		if(sum(null.perms)>0){
			xbeta.new <- lapply(1:perm.ind, function(j) lapply(1:k, function(i) X.aug.old[[i]] %*% matrix(beta.new[[j]][[i]],ncol=1)))
			res <- lapply(1:perm.ind, function(j) sapply(1:k, function(i) (y - xbeta.new[[j]][[i]])^2))
			s.new <- lapply(1:perm.ind, function(j) sqrt(sapply(1:k, function(i) sum(z[,perms[j,i]] * (res[[j]][, i]))/sum(z[,perms[j,i]]))))
			comp <- lapply(1:perm.ind, function(j) lapply(1:k, function(i) lambda.new[i] * dnorm(y, xbeta.new[[j]][[i]], s.new[[j]][i])))
			comp <- lapply(1:perm.ind, function(j) sapply(comp[[j]], cbind))
			compsum <- lapply(1:perm.ind, function(j) apply(comp[[j]], 1, sum))
			newobsloglik <- sapply(1:perm.ind, function(j) sum(log(compsum[[j]])))
			newobsloglik[c(1-null.perms)]=-Inf
			IND <- which.max(newobsloglik)
			z = z[,perms[IND,]]
			lambda.new <- apply(z, 2, mean)
			lambda <- lambda.new
			beta <- beta.new[[IND]]
			xbeta <- xbeta.new[[IND]]
			res <- res[[IND]]
			X.aug <- X.aug.old
			s <- s.new[[IND]]
			psi.locs <- psi.locs.old
			sing <- sum(s < 1e-08)
			newobsloglik <- newobsloglik[IND]
		} else{
			newobsloglik=Inf
			sing=1
		}
	}
	if(newobsloglik<obsloglik & !is.na(newobsloglik) & abs(newobsloglik) != Inf){
		z <- z.old
		lambda <- apply(z, 2, mean)
		X.aug.1 <- lapply(1:k, function(i) cbind(1,aug.x(x[,-1],unlist(psi.locs.old[[i]]),psi[i,],delta[[i]])))
		lm.out.old=lapply(1:k, function(i) lm(y~X.aug.1[[i]][,-1],weights=z[,i]))
		beta <- lapply(lm.out.old,coef)
		xbeta <- lapply(1:k, function(i) X.aug.1[[i]] %*% matrix(beta[[i]],ncol=1))
		res <- sapply(1:k, function(i) (y - xbeta[[i]])^2)
		s <- sqrt(sapply(1:k, function(i) sum(z[,i] * (res[, i]))/sum(z[,i])))
		psi.locs <- psi.locs.old
		comp <- lapply(1:k, function(i) lambda[i] * dnorm(y, xbeta[[i]], s[i]))
		comp <- sapply(comp, cbind)
		compsum <- apply(comp, 1, sum)
		newobsloglik <- sum(log(compsum))
		if(is.na(newobsloglik)){
			sing <- 1
		} else{
			if(newobsloglik<obsloglik) sing <- 1
		}
	}
	}
	if (sing > 0 || is.na(newobsloglik) || abs(newobsloglik) == Inf) {
		cat("Need new starting values due to singularity...", "\n")
		restarts <- restarts + 1
		if (restarts > max.restarts) stop("Too many tries!")
		tmp.ind=1
		while(tmp.ind==1){
			if(sum(1-sapply(delta,is.null))>0) psi.temp=psi.locs
			tmp <- try(suppressWarnings(segregmix.init(y=y, x=x.old, lambda = NULL, beta = NULL, s = NULL, k = k, seg.Z=seg.Z, 
				psi=psi, psi.locs = NULL)),silent=TRUE)
		if(class(tmp)!="try-error") tmp.ind=2
	}
	lambda <- tmp$lambda
	beta <- tmp$beta
	s <- tmp$s
	k <- tmp$k
	psi.locs <- tmp$psi.locs
	n <- length(y)
	diff <- 1
	iter <- 0
	X.aug <- lapply(1:k, function(i) cbind(1,aug.x(x[,-1],unlist(psi.locs[[i]]),psi[i,],delta[[i]])))
	xbeta <- lapply(1:k, function(i) X.aug[[i]] %*% matrix(beta[[i]],ncol=1))
	res <- sapply(1:k, function(i) as.numeric((y - xbeta[[i]])^2))
	comp <- t((lambda/sqrt(2 * pi * s^2)) * t(exp(-t(t(res)/(2 * s^2)))))
	obsloglik <- sum(log(apply(comp, 1, sum)))
	ll <- obsloglik
	} else {
		diff <- newobsloglik - obsloglik
		obsloglik <- newobsloglik
		ll <- c(ll, obsloglik)
		X.aug.old <- X.aug
		psi.locs.old <- psi.locs
		iter <- iter + 1
		if (verb) {
			cat("iteration =", iter, "diff =", diff, "log-likelihood =", obsloglik, "\n")
		}
	}
	}
	if (iter == maxit) {
		warning("Maximum number of iterations reached.", call. = FALSE)
	}
	if (iter == 1) {
		cat("Converged in 1 iteration. Consider rerunning with different starting values or smaller stopping criterion.", "\n")
	}
	cat("number of iterations=", iter, "\n")
	names(delta) <- c(paste("comp", ".", 1:k, sep = ""))
	names(seg.Z) <- c(paste("comp", ".", 1:k, sep = ""))
	names(psi.locs) <- c(paste("comp", ".", 1:k, sep = ""))
	names(beta) <- c(paste("comp", ".", 1:k, sep = ""))
	for(i in 1:k){
		names(beta[[i]])[1]="(Intercept)"
		names(beta[[i]])[2:ncol(x)]=colnames(x)[-1]
		if(!is.null(psi.locs[[i]])){
			for(j in 1:ncol(psi)){
				if(psi[i,j]>0 & j==1){
					names(beta[[i]])[(ncol(x)+1):(ncol(x)+cumsum(psi[i,])[j])]=c(paste(colnames(x)[j+1], ".", 1:psi[i,j], sep = ""))
				} else if(psi[i,j]>0) names(beta[[i]])[(ncol(x)+cumsum(psi[i,])[j-1]+1):(ncol(x)+cumsum(psi[i,])[j])]=c(paste(colnames(x)[j+1], ".", 1:psi[i,j], sep = ""))
			}
		}
	}
	colnames(z) <- c(paste("comp", ".", 1:k, sep = ""))
		a = list(x = x, y = y, lambda = lambda, beta = beta, 
		sigma = s, seg.Z = seg.Z, psi.locs = psi.locs, delta = delta, loglik = obsloglik, posterior = z, all.loglik = ll, 
		restarts = restarts, ft = "segregmixEM")
	class(a) = "mixEM"
	a
}