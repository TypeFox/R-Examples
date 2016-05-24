# Copyright (C) 2011 Pierrick Bruneau, see README for full notice



randomGmm <- function(domain=10) {
	# generates random 2D GMMs
	mod <- newGmm()
	K <- rpois(1, 5)
	if(K == 0) K <- rpois(1, 5)
	# generate weights from dirichlet
	mod$w <- rDirichlet(K)
	
	# generate means & covs from runif and covgen over some specified domain
	for(i in 1:K) {
		mod$mean[[i]] <- numeric()
		for(j in 1:2) {
			mod$mean[[i]][j] <- runif(1, min=-domain, max=domain)
		}
		
		mod$cov[[i]] <- covgen(bounds=c(2, 10))
	}
	return(mod)
}

plotGmm <- function(mod, steps=200) {
	# get domain from means and covariances
	K <- length(mod$w)
	if(length(mod$mean[[1]]) != 2) stop("only 2D GMMs can be used")
	
	low_left <- c(+Inf, +Inf)
	up_right <- c(-Inf, -Inf)
	
	for(i in 1:K) {
		cur_mean <- mod$mean[[i]]
		cur_vars <- diag(mod$cov[[i]])
		for(d in 1:2) {
			if((cur_mean[d] - 2 * cur_vars[d]) < low_left[d]) low_left[d] <- cur_mean[d] - 2 * cur_vars[d]
			if((cur_mean[d] + 2 * cur_vars[d]) > up_right[d]) up_right[d] <- cur_mean[d] + 2 * cur_vars[d]
		}
	}
	
	stepx <- (up_right[1] - low_left[1]) / steps
	stepy <- (up_right[2] - low_left[2]) / steps
	
	
	pdf <- matrix(nrow=steps+1, ncol=steps+1)
	
	for(i in 0:(steps+1)) {
		for(j in 0:(steps+1)) {
			pdf[i,j] <- gmmdensity(mod, t(as.matrix(c(low_left[1] + i * stepx, low_left[2] + j * stepy))))
		}
	}
	
	lattice::wireframe(pdf, shade = TRUE, 
		aspect = c(1, 1), 
		xlab="", ylab="", zlab="", light.source = c(10,0,10))


}


datagen <- function(dreal=2, deff=6, npts=200, noise=0.1, genmean=rep(0,dreal), genspan=6, iso=FALSE) {
	data <- matrix(0, nrow=npts, ncol=deff)
	
	combs <- rnorm(dreal * (deff - dreal), sd=1)
	combs <- matrix(combs, nrow=deff-dreal, ncol=dreal)
	
	gencov <- matrix(0, nrow=dreal, ncol=dreal)
	# random choice of diagonal elements
	for(i in 1:dreal) {
		gencov[i,i] <- runif(1, min=1, max=genspan)
	}

	noise <- genspan / 60

	# generate random corelations if isotropic is set to F
	if(!iso) {
		for(i in 1:dreal) {
			for(j in i:dreal) {
				if(i != j) {
					# through property found in wikipedia:Positive definiteness
					thres <- (gencov[i,i] + gencov[j,j]) / 2
					gencov[i,j] = runif(1, min=-thres, max=thres)
					gencov[j,i] = gencov[i,j]
				}
			}
		}
	}
	
	print(gencov)
	
	data[,1:dreal] <- mnormt::rmnorm(npts, mean=genmean, varcov=gencov)
	
	if(dreal < deff) {
		for(i in (dreal+1):deff) {
			data[,i] <- rnorm(npts, sd=noise)
			for(j in 1:dreal) {
				data[,i] <- data[,i] + combs[i-dreal, j] * data[,j]
			}
		}
	}
	
	return(list(data=data, genmean=genmean, gencov=gencov))
}





covgen <- function(d=2, bounds=c(1, 5)) {
	gencov <- matrix(0,nrow=d, ncol=d)
	done <- FALSE
	while(!done) {
		for(i in 1:d) {
			gencov[i,i] <- runif(1, min=bounds[1], max=bounds[2])
		}
	
		for(i in 1:d) {
			for(j in i:d) {
				if(i != j) {
					thres <- (gencov[i,i] + gencov[j,j]) / 2
					gencov[i,j] <- gencov[j,i] <- runif(1, min=-thres, max=thres)
				}
			}
		}
		# test for positive definiteness
		vals <- eigen(gencov, only.values=TRUE)$values
		cond <- vals[length(vals)] / vals[1]
		if(cond > 10^(-10)) done <- TRUE
	}
	return(gencov)
}




circlegen <- function(npts=200, radius=10, noise=1) {
	dat <- matrix(0, nrow=npts, ncol=2)
	
	for(i in 1:npts) {
		cur <- runif(1, min=0, max=2*pi)
		dat[i,] <- c(radius*cos(cur), radius*sin(cur)) + mnormt::rmnorm(1, varcov=noise*diag(2))
	}
	
	return(dat)
}

l2norm <- function(vec) {
	val <- 0
	for(i in 1:length(vec)) {
		val <- val + vec[i]^2
	}
	return(sqrt(val))
}

mppcaToGmm <- function(model, notau=FALSE) {
	# convert W nu... format to w mean cov
	output <- list()
	output$w <- numeric()
	output$mean <- list()
	output$cov <- list()
	
	d <- length(as.numeric(model$mumean[[1]]))
	
	output$w <- model$alpha
	output$w <- output$w / sum(output$w)
	for(i in 1:length(model$alpha)) {
		output$mean[[i]] <- model$mumean[[i]]
		if(!notau) {
			output$cov[[i]] <- model$wmean[[i]] %*% t(model$wmean[[i]]) + (1/model$taumoment[i]) * diag(d)
		} else {
			output$cov[[i]] <- model$wmean[[i]] %*% t(model$wmean[[i]]) + (10^(-1)) * diag(d)
		}
	}

	return(output)
}






gaussianKL <- function(N0, N1) {
	d <- length(N0[1,])
	return(log(det(N1)) - log(det(N0)) + sum(diag(solve(N1) %*% N0)) - d)
}

getQforComp <- function(loadings, tau=1.0, verbose=FALSE, quick=FALSE) {
	d <- length(loadings[,1])
	q <- length(loadings[1,])
	
	if(!quick) {
		cur <- loadings[,1] %*% t(loadings[,1]) + (1/tau) * diag(d)
		curind <- 1
		delta <- +Inf
		while((delta > 0.1) && (curind < q)) {
			curind <- curind + 1
			new <- loadings[,1:curind] %*% t(loadings[,1:curind]) + (1/tau) * diag(d)
			delta <- gaussianKL(new, cur)
			cur <- new

		}
		if(delta > 0.1) {
			# last dimension added still has an interest : curind is q.
			curind <- q
			if(verbose) print(paste("all dimensions are used, q is", curind))
		} else {
			if(verbose) print(paste("adding dim.", curind, "makes delta=", delta, ", q set to", curind-1))
			curind <- curind - 1
		}
	} else {
		# quicker solution : analyze column norms
		norms <- numeric()
		for(i in 1:q) {
			norms[i] <- sqrt(sum(loadings[,i] * loadings[,i]))
		}
		print(paste("min:", min(norms), "max:", max(norms)))
		tot <- sum(norms)
		cur <- 0
		curind <- 0
		while((cur/tot)<0.95) {
			curind <- curind + 1
			cur <- cur + norms[curind]
		}
	}
	return(curind)
	
}

isNonVoid <- function(loadings) {
	d <- length(loadings[,1])
	q <- length(loadings[1,])
	result <- FALSE
	for(i in 1:q) {
		for(j in 1:d) {
			if(loadings[j,i] > 0.1) result <- TRUE
		}
	}
	return(result)	
}


subVarbayes <- function(model, thres=2.001) {
	# "sub" a model returned by varbayes
	d <- length(model$data[1,])
	
	filter1 <- list()
	filter1$model <- list()
	filter1$model$alpha <- numeric()
	filter1$model$beta <- numeric()
	filter1$model$nu <- numeric()
	filter1$model$mean <- list()
	filter1$model$wish <- list()
	#filter1$model$resp <- numeric()
	filter1$data <- model$data
	
	selection <- numeric()
	for(i in 1:length(model$model$alpha)) {
		if(model$model$alpha[i] > thres) {
			selection <- c(selection, i)
			filter1$model$alpha <- c(filter1$model$alpha, model$model$alpha[i])
			filter1$model$nu <- c(filter1$model$nu, model$model$nu[i])
			filter1$model$beta <- c(filter1$model$beta, model$model$beta[i])
			filter1$model$mean <- appendToList(filter1$model$mean, model$model$mean[[i]])
			filter1$model$wish <- appendToList(filter1$model$wish, model$model$wish[[i]])
		}
	}
	filter1$model$resp <- model$model$resp[,selection]	
	return(filter1)
}


appendToMppca <- function(mppca1, mppca2) {
	# append mppca2 to mppca1
	# uses appendToList
	mppca1$alpha <- c(mppca1$alpha, mppca2$alpha)
	mppca1$numoment <- appendToList(mppca1$numoment, mppca2$numoment, appendList=TRUE)
	mppca1$nua <- c(mppca1$nua, mppca2$nua)
	mppca1$nub <- appendToList(mppca1$nub, mppca2$nub, appendList=TRUE)
	mppca1$taumoment <- c(mppca1$taumoment, mppca2$taumoment)
	mppca1$taua <- c(mppca1$taua, mppca2$taua)
	mppca1$taub <- c(mppca1$taub, mppca2$taub)
	mppca1$wmean <- appendToList(mppca1$wmean, mppca2$wmean, appendList=TRUE)
	mppca1$wsigma <- appendToList(mppca1$wsigma, mppca2$wsigma, appendList=TRUE)
	mppca1$xsigma <- appendToList(mppca1$xsigma, mppca2$xsigma, appendList=TRUE)
	mppca1$mumean <- appendToList(mppca1$mumean, mppca2$mumean, appendList=TRUE)
	mppca1$musigma <- appendToList(mppca1$musigma, mppca2$musigma, appendList=TRUE)
	mppca1$mustar <- appendToList(mppca1$mustar, mppca2$mustar, appendList=TRUE)
	
	return(mppca1)
}




subMppca <- function(model, prune=FALSE, thres=2.001, quick=FALSE, noxmean=TRUE) {
	# ncol for selecting default number of column factors
	
	d <- length(as.numeric(model$mumean[[1]]))
	q <- length(model$numoment[[1]])

	# select xmean, x1mean or x2mean, depending on cases.

	filter1 <- list()
	filter1$alpha <- numeric()
	filter1$numoment <- list()
	filter1$nua <- numeric()
	filter1$nub <- list()
	filter1$taumoment <- numeric()
	filter1$taua <- numeric()
	filter1$taub <- numeric()
	filter1$wmean <- list()
	filter1$wsigma <- list()
	filter1$xsigma <- list()
	filter1$mumean <- list()
	filter1$musigma <- list()
	filter1$mustar <- list()
	if(!is.null(model$xmean)) {
		if(!noxmean) filter1$xmean <- list()
	} else if(!is.null(model$x1mean)) {
		if(!noxmean) {
			filter1$x1mean <- list()
			filter1$x2mean <- list()
		}
	}



	for(i in 1:length(model$alpha)) {
		if(model$alpha[i] > thres) {
			filter1$alpha <- c(filter1$alpha, model$alpha[i])
			filter1$numoment <- appendToList(filter1$numoment, model$numoment[[i]])
			filter1$nua <- c(filter1$nua, model$nua[i])
			filter1$nub <- appendToList(filter1$nub, model$nub[[i]])
			filter1$taumoment <- c(filter1$taumoment, model$taumoment[i])
			filter1$taua <- c(filter1$taua, model$taua[i])
			filter1$taub <- c(filter1$taub, model$taub[i])
			filter1$wmean <- appendToList(filter1$wmean, model$wmean[[i]])
			filter1$wsigma <- appendToList(filter1$wsigma, model$wsigma[[i]])
			filter1$xsigma <- appendToList(filter1$xsigma, model$xsigma[[i]])
			filter1$mumean <- appendToList(filter1$mumean, model$mumean[[i]])
			filter1$musigma <- appendToList(filter1$musigma, model$musigma[[i]])
			filter1$mustar <- appendToList(filter1$mustar, model$mustar[[i]])
			if(!is.null(filter1$xmean)) {
				filter1$xmean <- appendToList(filter1$xmean, model$xmean[[i]])
			} else if(!is.null(filter1$x1mean)) {
				filter1$x1mean <- appendToList(filter1$x1mean, model$x1mean[[i]])
				filter1$x2mean <- appendToList(filter1$x2mean, model$x2mean[[i]])
			}

		}
	}
	


	qmax <- 1
	

	
	if(prune) {
		for(i in 1:length(filter1$alpha)) {
			if(isNonVoid(filter1$wmean[[i]])) {
				curq <- getQforComp(filter1$wmean[[i]], filter1$taumoment[i], quick=quick)
				if(curq > qmax) qmax <- curq 
			}
		}
	} else {
		qmax <- q
	}
	
	print(paste(qmax, "factors selected"))

	

	filter2 <- list()

	filter2$alpha <- filter1$alpha
	filter2$numoment <- list()
	filter2$nua <- filter2$nua
	filter2$nub <- list()
	filter2$taumoment <- filter1$taumoment
	filter2$taua <- filter1$taua
	filter2$taub <- filter1$taub
	filter2$wmean <- list()
	filter2$wsigma <- list()
	filter2$xsigma <- list()
	filter2$mumean <- filter1$mumean
	filter2$musigma <- filter1$musigma
	filter2$mustar <- filter1$mustar
	if(!is.null(filter1$xmean)) {
		if(!noxmean) filter2$xmean <- list()
	} else if(!is.null(filter1$x1mean)) {
		if(!noxmean) {
			filter2$x1mean <- list()
			filter2$x2mean <- list()
		}
	}


	for(i in 1:length(filter1$alpha)) {
		filter2$numoment[[i]] <- filter1$numoment[[i]][1:qmax]
		filter2$nub[[i]] <- filter1$nub[[i]][1:qmax]
		if(qmax > 1) {
			filter2$wmean[[i]] <- filter1$wmean[[i]][,1:qmax]
			filter2$wsigma[[i]] <- filter1$wsigma[[i]][1:qmax, 1:qmax]
			filter2$xsigma[[i]] <- filter1$xsigma[[i]][1:qmax, 1:qmax]
			if(!is.null(filter2$xmean)) {
				filter2$xmean[[i]] <- filter1$xmean[[i]][,1:qmax]
			} else if(!is.null(filter2$x1mean)) {
				filter2$x1mean[[i]] <- filter1$x1mean[[i]][,1:qmax]
				filter2$x2mean[[i]] <- list()
				for(j in 1:length(filter1$x2mean[[i]])) {
					filter2$x2mean[[i]][[j]] <- filter1$x2mean[[i]][[j]][1:qmax, 1:qmax]
				}
			}
		} else {
			filter2$wmean[[i]] <- as.matrix(filter1$wmean[[i]][,1])
			filter2$wsigma[[i]] <- filter1$wsigma[[i]][1,1]
			filter2$xsigma[[i]] <- filter1$xsigma[[i]][1,1]
			if(!is.null(filter2$xmean)) {
				filter2$xmean[[i]] <- as.matrix(filter1$xmean[[i]][,1])
			} else if(!is.null(filter2$x1mean)) {
				filter2$x1mean[[i]] <- as.matrix(filter1$x1mean[[i]][,1])
				filter2$x2mean[[i]] <- list()
				for(j in 1:length(filter1$x2mean[[i]])) {
					filter2$x2mean[[i]][[j]] <- filter1$x2mean[[i]][[j]][1,1]
				}

			}	

		}
	}

	return(filter2)
	
}



incremMerge <- function(modref, newmod, k=200, nit=100, quick=FALSE) {
	modref$alpha <- modref$alpha * 100000 / sum(modref$alpha)
	newmod$alpha <- newmod$alpha * 100000 / sum(newmod$alpha)
	# inputs are assumed to be correctly filtered
	
	modref <- appendToMppca(modref, newmod)
	modref <- normMppca(modref)

	modref <- mmppca(modref, k, maxit=nit)
	modref <- subMppca(modref)

	return(modref)
}



ZtoLabels <- function(resp) {
	labs <- numeric()
	for(i in 1:length(resp[,1])) {
		labs[i] <- which.max(resp[i,])
	}
	return(labs)
}





subGmm <- function(model, dims=c(1,2), inds=NULL) {
	output <- list()
	output$mean <- list()
	output$cov <- list()

	
	if(is.null(inds)) inds <- 1:length(model$w)
	output$w <- model$w[inds]
	output$w <- output$w / sum(output$w)
	
	#for(i in 1:length(model$w)) {
	for(i in 1:length(inds)) {
		cur <- model$mean[[inds[i]]]
		cur <- as.numeric(cur)
		cur <- cur[dims]
		cur <- as.matrix(cur)
		output$mean[[i]] <- cur
		output$cov[[i]] <- model$cov[[inds[i]]][dims, dims]
	}
	return(output)
}


newMppca <- function() {
	out <- list()
	out$alpha <- numeric()
	out$wmean <- list()
	out$mumean <- list()
	out$numoment <- list()
	out$taumoment <- numeric()
	return(out)
}


gmmToMppca <- function(model, alpha=500) {
	# alpha is the overall population assicated with output MPPCA
	output <- newMppca()
	k <- length(model$w)
	d <- length(as.numeric(model$mean[[1]]))
	
	q <- 1
	eigens <- list()
	# get the dimensionality to retain
	for(i in 1:k) {
		eigens[[i]] <- eigen(model$cov[[i]])
		ratio <- 0
		curq <- 0
		while(ratio < 0.8) {
			curq <- curq + 1
			ratio <- sum(eigens[[i]]$values[1:curq]) / sum(eigens[[i]]$values)
		}
		if(curq > q) q <- curq		
	}
	
	print(paste("q selected :", q))
	
	for(i in 1:k) {
		output$alpha[i] <- model$w[i] * alpha
		output$mumean[[i]] <- model$mean[[i]]
		if(q > 1) {
			output$wmean[[i]] <- eigens[[i]]$vectors[,1:q] %*% diag(eigens[[i]]$values[1:q])
		} else {
			output$wmean[[i]] <- as.matrix(eigens[[i]]$vectors[,1]) * eigens[[i]]$values[1]
		}
		# set columns to same sign
		for(j in 1:q) {
			if(output$wmean[[i]][1,j] < 0) output$wmean[[i]][,j] <- - output$wmean[[i]][,j]
		}
		
		output$numoment[[i]] <- 1 / eigens[[i]]$values[1:q]
		output$taumoment[i] <- 1 / mean(eigens[[i]]$values[(q+1):d])
	}
	
	return(output)
}



newGmm <- function() {
	newmod <- list()
	newmod$w <- numeric()
	newmod$mean <- list()
	newmod$cov <- list()
	newmod$a <- numeric()
	return(newmod)
}

appendToGmm <- function(mod1, mod2) {
	if(length(mod1$a)==0) {
		mod2$a <- rep(0, length(mod2$w))
	} else {
		mod2$a <- rep(mod1$a[length(mod1$a)] + 1, length(mod2$w))
	}
	mod1$a <- c(mod1$a, mod2$a)
	mod1$w <- c(mod1$w, mod2$w)
	mod1$mean <- appendToList(mod1$mean, mod2$mean, appendList=TRUE)
	mod1$cov <- appendToList(mod1$cov, mod2$cov, appendList=TRUE)
	return(mod1)
}



normMppca <- function(mppca1) {
	maxq <- 0
	len <- length(mppca1$alpha)
	for(i in 1:len) {
		curq <- length(mppca1$numoment[[i]])
		if(curq > maxq) maxq <- curq
	}
	for(i in 1:len) {
		curq <- length(mppca1$numoment[[i]])
		d <- length(mppca1$wmean[[i]][,1])
		if(curq < maxq) {
			newwmean <- matrix(0, nrow=d, ncol=maxq)
			newwmean[,1:curq] <- mppca1$wmean[[i]]
			mppca1$wmean[[i]] <- newwmean
			mppca1$numoment[[i]][(curq+1):maxq] <- 100
		}
	}
	print(paste("normalized at", maxq, "factors"))
	return(mppca1)
}



eigenMppca <- function(mod) {
	# diagonalize t(W) %*% W => eigenvectors (cols) are t(R)
	# input is assumed to having been reduced (subMppca)
	
	d <- length(mod$mumean[[1]])
	
	for(i in 1:length(mod$alpha)) {
		# first build the matrix to diagonalize
		mat <- t(mod$wmean[[i]]) %*% mod$wmean[[i]]
		
		# perform the diag
		op <- eigen(mat)
		Rt <- op$vectors
		vals <- op$values
		
		# see tipping PPCA
		mod$wmean[[i]] <- mod$wmean[[i]] %*% Rt
		
	}
	return(mod)
	
}

gramschmidt <- function(mat) {
	# perform gramschmidt normalisation on the set of column vectors in mat
	
	d <- length(mat[1,])
	for(i in 1:d) {
		for(j in 1:i) {
			if(j < i) {
				# previous vectors are normalized : no norm term.
				proj <- mat[,j] * as.numeric(t(mat[,i]) %*% mat[,j])
				mat[,i] <- mat[,i] - proj
			}
		}
		# when finished, normalize current vector
		norm <- as.numeric(t(mat[,i]) %*% mat[,i])
		mat[,i] <- mat[,i] / sqrt(norm)
	}
	return(mat)
}


semispheregen <- function(npts=200, radius=10, noise=1) {
	dat <- matrix(nrow=npts, ncol=3)
	
	for(i in 1:npts) {
		theta <- runif(1, min=0, max=pi/2)
		phi <- runif(1, min=-pi, max=pi)
		dat[i,] <- c(radius*cos(theta)*cos(phi), radius*cos(theta)*sin(phi), radius*sin(theta)) + 
			mnormt::rmnorm(1, varcov=noise*diag(3))
	}
	return(dat)
}







normalizeVariable <- function(v) {
	themin <- min(v)
	themax <- max(v)
	v <- v - themin
	v <- v / (themax - themin)
	return(v)
}


binnedEntropy <- function(v, nbins=100) {
	v <- normalizeVariable(v)
	step <- 1/nbins
	counts <- numeric()
	cumul <- 0
	
	for(i in 1:nbins) {
		cur <- sum(as.numeric(v <= step*i)) - cumul
		counts <- c(counts, cur)
		cumul <- cumul + cur
	}
	counts <- counts / sum(counts)
	
	entropy <- 0
	for(i in 1:nbins) {
		if(counts[i] != 0) {
			entropy <- entropy - counts[i] * log(counts[i]) 
		}
	}
	return(entropy)
}


reBuild <- function(v, voids, nonvoids, domains, placeholder=1) {
	# rescale v to correct domain
	v <- t(as.matrix(as.numeric(v)))
	
	v <- as.numeric(setDomain(v, domains, 10))
	v2 <- numeric()
	for(i in 1:length(v)) {
		v2[nonvoids[i]] <- v[i]
	}
	
	v2[voids] <- placeholder
	v2 <- pixmap::pixmapGrey(matrix(v2, nrow=sqrt(length(v2)), ncol=sqrt(length(v2))))

	return(v2)
	
}



pixmapToVector <- function(p) {
	return(as.numeric(pixmap::getChannels(p, "grey")))
}



spiralgen <- function(radius=10, n=1000, laps=2, noise=1) {
	nptsperlap <- n / laps
	res <- matrix(nrow=n, ncol=2)
	for(i in 1:n) {
		curangle <- ((i %% nptsperlap) / nptsperlap) * (2 * pi)
		curradius <- (1+ (i / nptsperlap)) * (radius)
		res[i,1] <- curradius * cos(curangle) + rnorm(1, sd=noise)
		res[i,2] <- curradius * sin(curangle) + rnorm(1, sd=noise)
	}
	return(res)
	
}








setDomain <- function(dat, span=10, oldspan=NULL) {
	# possible lengths for span :
	# - class = numeric et 1 = all axes on [-span, span]
	# - class = numeric and 2 = all axes on [span[1], span[2]]
	# - class = list et length 2*d = full specification of axes.
	# see "set domain demonstration"
	
	n <- length(dat[,1])
	d <- length(dat[1,])
	
	maxs2 <- numeric()
	mins2 <- numeric()
	maxs1 <- numeric()
	mins1 <- numeric()
	
	# symmetric or asymetric domains ; depends on length of "span" argument
	if((class(span) == "numeric") && (length(span) == 1)) {
		if(span < 0) span <- -span
		maxs2 <- rep(span, d)
		mins2 <- rep(-span, d)
	} else if(class(span) == "numeric" && (length(span) == 2)) {
		if(span[1] > span[2]) {
			temp <- span[1]
			span[1] <- span[2]
			span[2] <- temp
		}
		maxs2 <- rep(span[2], d)
		mins2 <- rep(span[1], d)
	} else if((class(span) == "list") && (length(span) == 2)) {
		for(i in 1:d) {
			if(span[[1]][i] > span[[2]][i]) {
				temp <- span[[1]][i]
				span[[1]][i] <- span[[2]][i]
				span[[2]][i] <- temp
			}
			maxs2[i] <- span[[2]][i]
			mins2[i] <- span[[1]][i]
		}
	} else {
		stop("inadequate 'span' argument")
	}
	
	

	if(is.null(oldspan)) {
		mins1 <- apply(dat, 2, min)
		maxs1 <- apply(dat, 2, max)		
	} else if((class(oldspan) == "numeric") && (length(oldspan) == 1)) {
		if(oldspan < 0) oldspan <- -oldspan
		maxs1 <- rep(oldspan, d)
		mins1 <- rep(-oldspan, d)
	} else if(class(oldspan) == "numeric" && (length(oldspan) == 2)) {
		if(oldspan[1] > oldspan[2]) {
			temp <- oldspan[1]
			oldspan[1] <- oldspan[2]
			oldspan[2] <- temp
		}
		maxs1 <- rep(oldspan[2], d)
		mins1 <- rep(oldspan[1], d)
	} else if((class(oldspan) == "list") && (length(oldspan) == 2)) {
		for(i in 1:d) {
			if(oldspan[[1]][i] > oldspan[[2]][i]) {
				temp <- oldspan[[1]][i]
				oldspan[[1]][i] <- oldspan[[2]][i]
				oldspan[[2]][i] <- temp
			}
			maxs1[i] <- oldspan[[2]][i]
			mins1[i] <- oldspan[[1]][i]
		}
	} else {
		stop("inadequate 'oldspan' argument")
	}
	
	
	
	
	for(i in 1:n) {
		for(j in 1:d) {
			dat[i,j] <- (mins2[j] * (maxs1[j] - dat[i,j]) + maxs2[j] * (dat[i,j] - mins1[j])) / (maxs1[j] - mins1[j]) 
		}
	}
	
	return(dat)
}




pca <- function(dat, ncomp=NULL) {
	# force-cast to matrix
	dat <- as.matrix(dat)

	# if ncomp is NULL, all columns are selected
	if(is.null(ncomp)) ncomp <- dim(dat)[2]
	
	# nombre d'individus = longueur d'une colonne
	n <- length(dat[,1])
	
	# matrice des poids, vecteur de 1
	D <- diag(rep(1/n, n))
	one <- as.matrix(rep(1,n))
	
	# centred data
	Y <- (diag(n) - one %*% t(one) %*% D) %*% dat
	
	# variance
	V <- t(Y) %*% D %*% Y
	
	# D(1/s2) metric
	d1s2 <- diag(1/diag(V))
	
	# diagonalisation
	eig <- eigen(V %*% d1s2)

	eigvecs <- Re(eig$vectors)
	
	# fonction pour M-normer chaque vecteur de notre ensemble de vecteurs propres
	# en effet eigen les I-norme a 1, nous on veut utiliser D(1/s2)
	A <- apply(eigvecs, 2, function(X) { norm <- sqrt(as.numeric(t(X) %*% d1s2 %*% X))
	return(X / norm)})
	
	#info <- "eigenvalues :"
	#for(i in 1:length(eig$values)) {
	#	info <- paste(info, eig$values[i])
	#}
	#
	#print(info)
	eigvals <- Re(eig$values)
	message("ncomp=", ncomp, ", with ", format(100 * sum(eigvals[1:ncomp]) / sum(eigvals), digits=2),
		"% of retained variance")
	
	# U = transposee(A^{-1}), cf cours et TP ACP
	U <- t(solve(A))
	
	# transformed data (XU dans le TP)
	reduce <- Y %*% U
	
	# on retourne les ncomp premieres composantes (2D par defaut)
	return(reduce[,1:ncomp])
}



getDataLikelihood <- function(gmm, dat) {
	densities <- gmmdensity(gmm, dat)
	res <- sum(log(densities))
	return(res)
}


getBic <- function(gmm, dat) {
	lnL <- getDataLikelihood(gmm, dat)
	n <- length(dat[,1])
	d <- length(dat[1,])
	k <- length(gmm$w)
	# constraint on w
	paramsize <- (k-1) + k*d + k*(d*(d+1))/2
	res <- -2 * lnL + paramsize * log(n)
	return(res)
}


getVarbayesResp <- function(data, model) {
	# get responsibilities according to varbayes model
	n <- dim(data)[1]
	d <- dim(data)[2]
	k <- length(model$model$alpha)
	resp <- matrix(nrow=n, ncol=k)
	gammamat <- matrix(nrow=k, ncol=d)
	kappalog <- numeric()
	
	for(j in 1:d) {
		#print(model$model$nu)
		#print(digamma((model$model$nu - j + 1) / 2))
		#print(k)
		gammamat[,j] <- digamma((model$model$nu - j + 1) / 2)
	}
	
	for(j in 1:k) {
		kappalog[j] <- sum(gammamat[j,]) + 2 * d + log(det(model$model$wish[[j]]))
	}
	
	for(i in 1:n) {
		resp[i,] <- digamma(model$model$alpha) - digamma(sum(model$model$alpha))
		resp[i,] <- resp[i,] + 0.5 * kappalog
		for(j in 1:k) {
			resp[i,j] <- resp[i,j] - 0.5 * (d / model$model$beta[j] + model$model$nu[j] * t(as.matrix(data[i,] - model$model$mean[[j]])) %*% model$model$wish[[j]] %*% as.matrix(data[i,] - model$model$mean[[j]]))
		}
		themax <- max(resp[i,])
		resp[i,] <- resp[i,] - themax
		resp[i,] <- exp(resp[i,])
		resp[i,] <- resp[i,] / sum(resp[i,])
	}
	return(resp)
}


appendToList <- function(lst, obj, appendList=FALSE) {
	if(appendList) {
		for(i in 1:length(obj)) {
			lst[[length(lst) + 1]] <- obj[[i]]
		}
	} else {
		lst[[length(lst) + 1]] <- obj
	}
	return(lst)
}



generate2Dtransform <- function(dims=4) {
	transform <- matrix(rnorm(2*dims), nrow=dims, ncol=2)
	# check linear independency <=> X^T %*% X has full rank
	
	while(det(t(transform) %*% transform) < 0.001) {
		transform <- matrix(rnorm(2*dims), nrow=dims, ncol=2)
		message("degenerate 2D transform - trying again...")
	}
	
	transform <- gramschmidt(transform)
	# 2 first dims = original signal
	transform[1:2, 1:2] <- diag(2)
	return(transform)
}



dat1sample <- function(nelts, gmm, noise, transform=generate2Dtransform(2), oldbounds=NULL, newbounds=NULL) {
	# generate elements
	d <- dim(transform)[1]
	dat <- gmmgen(gmm, nelts)[[1]]
	dat <- dat %*% t(transform) + matrix(rnorm(nelts * d, sd=noise), nrow=nelts, ncol=d)
	
	# change domain if necessary
	if(!is.null(newbounds)) {
		if(is.null(oldbounds)) {
			oldbounds <- list()
			oldbounds[[1]] <- apply(dat, 2, min)
			oldbounds[[2]] <- apply(dat, 2, max)
		}
		dat <- setDomain(dat, newbounds, oldbounds)
	}
	return(dat)
}



dat2sample <- function(nelts, radius, noise, oldbounds=NULL, newbounds=NULL) {
	dat <- semispheregen(nelts, radius, noise)
	
	# change domain if necessary
	if(!is.null(newbounds)) {
		if(is.null(oldbounds)) {
			oldbounds <- list()
			oldbounds[[1]] <- apply(dat, 2, min)
			oldbounds[[2]] <- apply(dat, 2, max)
		}
		dat <- setDomain(dat, newbounds, oldbounds)
	}
	return(dat)
	
}

dat3sample <- function(nelts, radius, noise, transform=generate2Dtransform(2), oldbounds=NULL, newbounds=NULL) {
	# noise for transform = circle noise / 10
	d <- dim(transform)[1]
	dat <- circlegen(nelts, radius, noise)
	dat <- dat %*% t(transform) + matrix(rnorm(nelts * d, sd=noise/10), nrow=nelts, ncol=d)	
	
	# change domain if necessary
	if(!is.null(newbounds)) {
		if(is.null(oldbounds)) {
			oldbounds <- list()
			oldbounds[[1]] <- apply(dat, 2, min)
			oldbounds[[2]] <- apply(dat, 2, max)
		}
		dat <- setDomain(dat, newbounds, oldbounds)
	}
	return(dat)
	
}

generateSparsePoints <- function(npoints, dim=2, span=10, mindist=2, maxit=20) {
	converge <- FALSE
	thepoints <- matrix(nrow=npoints, ncol=dim)
	notfar <- 1:npoints
	if(span < 0) span <- -span
	nit <- 1
	while((!converge) && (nit < maxit)) {
		print(paste("it", nit, ",", length(notfar), "points to gen."))	
		thepoints[notfar,] <- runif(length(notfar) * dim, -span, span)
		notfar <- which(.Call("control", thepoints, mindist, PACKAGE="VBmix") == 0)
		converge <- (length(notfar) == 0)
		nit <- nit + 1
	}
	
	if(!converge) print("could not converge, current points are returned")
	return(thepoints)
}






