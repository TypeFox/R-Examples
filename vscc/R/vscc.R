####Fix the arguments passed along in the post-analysis!

vscc <-function(x, G=1:9, automate="mclust", initial=NULL, train=NULL, forcereduction=FALSE){
	origx <- x
	origG <- G
	x <- scale(x)
	p <- ncol(x)
	if(is.null(train)){
		if(automate=="teigen"){
#			require("teigen")
			if(packageVersion("teigen")<1.9){
				warning(paste("The 'vscc' package requires 'teigen' version 1.9 or higher, version", packageVersion("teigen"), "is currently installed: issues may arise."))
			}
			if(is.null(initial)){
				mclinit1 <- hc("VVV",x)
				mclinit2 <- hclass(mclinit1, G)
				mclist <- list()
				for(g in length(G)){
					mclist[[G[g]]] <- mclinit2[,g]
				}
				initrun <- teigen(x, G, init=mclist, training=train, verbose=FALSE) 
				if(initrun$G==1){stop("teigen initialization gives G=1 solution...please use an initialization where G>1")}
				initial <- initrun$classification
				initunc <- sum(1-apply(initrun$fuzzy,1,max))
			}
			G <- length(unique(initial))
			n <- nrow(x)
			zmat <- matrix(0,n,G)
			for(i in 1:G){
				zmat[initial==i, i]<-1
			}
		}
		else{
			if(automate=="mclust"){
#				require("mclust")
				if(packageVersion("mclust")<4.0){
					warning(paste("VSCC requires 'mclust' version 4.0 or higher, version", packageVersion("mclust"), "is currently installed: issues may arise."))
				}
				if(is.null(initial)){
					initrun <- Mclust(x, G)
					if(initrun$G==1){stop("mclust initialization gives G=1 solution...please use an initialization where G>1")}
					initial <- initrun$classification
					initunc <- sum(initrun$unc)
				}
				G <- length(unique(initial))
				n <- nrow(x)
				zmat <- matrix(0,n,G)
				for(i in 1:G){
					zmat[initial==i, i]<-1
				}
			}
			else{
				if(is.null(initial)){stop("If an initial clustering vector is not supplied, automate='teigen' or 'mclust' must be specified")}
			}
			G <- length(unique(initial))
			n <- nrow(x)
			zmat <- matrix(0,n,G)
			for(i in 1:G){
				zmat[initial==i, i]<-1
			}
		}
	}
	else{
		if(is.null(initial)){
			stop("If using 'train', 'initial' vector must also be given")
		}
		origx <- x
		x <- x[train,]
		G <- length(unique(initial[train]))
		n <- nrow(x)
		zmat <- matrix(0,n,G)
		for(i in 1:G){
			zmat[initial[train]==i, i]<-1
		}
	}
	
	if(is.null(colnames(x))){ 
		colnames(x) <- 1:p 
		colnames(origx) <- 1:p
	}
	ng <- colSums(zmat)
	mug <- matrix(0,G,p)
	for(g in 1:G){
		mug[g,] <- colSums(zmat[,g]*x)/ng[g]
	}
#	mug <- muginit(G,p,x,zmat,ng)	
	
	mugarr <- array(0,dim=c(n,p,G))
	for(g in 1:G){
		mugarr[,,g] <- t(mug[g,] * t(matrix(1,n,p)))
	}
	mugmat <- matrix(0,n,p)
	for(g in 1:G){
		mugmat <- mugmat + zmat[,g] * mugarr[,,g]
	}
	xminusmug <- x - mugmat
	ss <- xminusmug * xminusmug
	ssbyvar <- colSums(ss)/n
	
#	bssmugmat <- array(0,dim=c(n,p,G))
#	for(g in 1:G){
#		bssmugmat[,,g] <- bssmugmat[,,g] + (1-zmat[,g]) * mugarr[,,g]
#	}
#	bssxminusmug <- array(0,dim=c(n,p,G))
#	for(g in 1:G){
#		bssxminusmug[,,g] <- (x-bssmugmat[,,g])^2/(n-ng[g])
#	}
#	bssbyvar <- rep(0,p)
#	for(g in 1:G){
#		bssbyvar <- bssbyvar + colSums(bssxminusmug[,,g])
#	}
#	sortbss <- sort(bssbyvar)
	sorted <- t(as.matrix(sort(ssbyvar)))
	select <- list()
	useselect <- list()
	varnames <- list()
	trun <- list()
	numvars <- NA
	for(i in 1:5){
		select[[i]] <-  matrix(data=origx[,colnames(sorted)[1]])
		useselect[[i]] <- matrix(data=x[,colnames(sorted)[1]])
		varnames[[i]] <- colnames(sorted)[1]
	}
	counts <- rep(2,5)
	for(k in 2:p){
		curname <- colnames(sorted)[k]
		for(i in 1:5){
			curcor <- cor(cbind(x[,curname],useselect[[i]]))
			if(all(abs(curcor[upper.tri(curcor)])<=(1-sorted[1,k]^i))){
				select[[i]] <- cbind(select[[i]],origx[,curname])
				useselect[[i]] <- cbind(useselect[[i]],x[,curname])
				varnames[[i]][counts[i]] <- curname
				counts[i] <- counts[i]+1
			}
		}
	}
	for(i in 1:5){
		colnames(select[[i]]) <- varnames[[i]]
	}
	if(!is.null(automate)){
		tuncs <- Inf
		numvars <- counts-1
		counttab <- table(counts-1)
		runteig <- rep(TRUE,5)
		if(any(counttab>1)){
#			dubvars <- as.numeric(names(which(counttab>1)))
#			for(j in 1:length(dubvars)){
#				relneedcheck <- which(numvars==dubvars[j])
#				k <- 1
#				while(k < length(relneedcheck)){
#					for(i in (k+1):length(relneedcheck)){
#						if(all(varnames[[relneedcheck[k]]] %in% varnames[[relneedcheck[i]]])){
#							runteig[relneedcheck[i]] <- FALSE
#						}
#					}
#					k <- k+1
#				}
#			}
			#This could be improved
			for(i in 1:4){
				for(j in (i+1):5){
					if(length(varnames[[i]])==length(varnames[[j]])){
						if(all(varnames[[i]] %in% varnames[[j]])){
							runteig[j] <- FALSE
						}
					}
				}
			}
		}
		if(automate=="teigen"){
			for(i in 1:5){
				if(runteig[i]){
					G <- origG
					mclinit1 <- hc("VVV",x)
					mclinit2 <- hclass(mclinit1, G)
					mclist <- list()
					for(g in length(G)){
						mclist[[G[g]]] <- mclinit2[,g]
					}
					trun[[i]] <- teigen(select[[i]], G, training=train, init=mclist, verbose=FALSE)
					if(trun[[i]]$G>1){
						tuncs[i] <- sum(1-apply(trun[[i]]$fuzzy,1,max))
					}
					else{
						tuncs[i] <- Inf
					}
				}
				else{
					trun[[i]] <- "Same as simpler relation"
					tuncs[i] <- Inf
				}
			}
		}
		else{
			if(is.null(train)){
				for(i in 1:5){
					if(runteig[i]){
						G <- origG
						trun[[i]] <- Mclust(scale(select[[i]]), G)
						if(trun[[i]]$G>1){
							tuncs[i] <- sum(trun[[i]]$unc)
						}
						else{
							tuncs[i] <- Inf
						}
					}
					else{
						trun[[i]] <- "Same as simpler relation"
						tuncs[i] <- Inf
					}
				}
			}
			else{
				for(i in 1:5){
					if(runteig[i]){
						G <- origG
						trun[[i]] <- teigen(select[[i]], G, models="mclust", training=train, init="uniform", verbose=FALSE, known=initial)
						if(trun[[i]]$G>1){
							tuncs[i] <- sum(1-apply(trun[[i]]$fuzzy,1,max))
						}
						else{
							tuncs[i] <- Inf
						}
					}
					else{
						trun[[i]] <- "Same as simpler relation"
						tuncs[i] <- Inf
					}
				}
			}
		}
		
	}
	store <- list()
	store[["selected"]] <- select
	if(!is.null(automate)){
#		if(is.null(initial)){
			if(is.null(train)){
				store[["initialrun"]] <- initrun
			}
			else{
				initunc <- Inf
			}
			if(forcereduction){
				store[["bestmodel"]] <- trun[[which.min(tuncs)]]
				store[["chosenrelation"]] <- which.min(tuncs)
			}
			else{
				if(min(tuncs)<initunc){
					store[["bestmodel"]] <- trun[[which.min(tuncs)]]
					store[["chosenrelation"]] <- which.min(tuncs)
					store[["topselected"]] <- select[[which.min(tuncs)]]
					store[["uncertainty"]] <- min(tuncs)
				}
				else{
					store[["bestmodel"]] <- initrun
					store[["chosenrelation"]] <- "Full dataset"
					store[["topselected"]] <- origx
					store[["uncertainty"]] <- initunc
				}
			}
#		}
#		else{
#			store[["bestmodel"]] <- trun[[which.min(tuncs)]]
#			store[["chosenrelation"]] <- which.min(tuncs)
#			store[["topselected"]] <- select[[which.min(tuncs)]]
#			store[["uncertainty"]] <- min(tuncs)
#		}
		store[["allmodelfit"]] <- trun
	}
	store[["family"]] <- automate
	store[["wss"]] <- sorted
	class(store) <- "vscc"
	store
}
