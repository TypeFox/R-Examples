
#prcomp$sdev = svd$d / sqrt(n-1)
MetaPCA <- function(DList, method=c("Angle","Eigen","RobustAngle","SparseAngle"), robust.var=c("qn","mad"), nPC=2,
		.weight=rep(1/length(DList),length(DList)), sparse.maxFeatures=NULL, sparse.lambda=NULL, sparse.max.iter=100, sparse.eps=1e-3,
		.scale=FALSE, .scaleAdjust=TRUE, doPreprocess=TRUE, cutRatioByMean=.4, cutRatioByVar=.4, doImpute=TRUE,
		na.rm.pct=.1, na.rm.pct.each=.5, verbose=FALSE) 
{
	method <- match.arg(method)
	stopifnot(nPC < min(sapply(DList, ncol)))
	
	if(method=="SparseAngle") {
		.res <- .MetaSparsePCA(DList, nPC=nPC, maxFeatures=sparse.maxFeatures, lambda=sparse.lambda, .scale=.scale, 
				.scaleAdjust=.scaleAdjust, max.iter=sparse.max.iter, eps=sparse.eps)
	} else {
		robust.var <- match.arg(robust.var)
		
		if(doPreprocess)
			DList <- PreprocessMetaAnalysis(DList, cutRatioByMean=cutRatioByMean, cutRatioByVar=cutRatioByVar, doImpute=doImpute,
					na.rm.pct=na.rm.pct, na.rm.pct.each=na.rm.pct.each, verbose=FALSE)
		
		for(i in 1:length(DList)) {
			DList[[i]] <- t(scale(t(DList[[i]]), scale=.scale))
		}
		
		.weight <- .weight / sum(.weight)
		
		.totalVar <- NULL
		.CovSum <- 0
		for(i in 1:length(DList)) {
			if(method=="RobustAngle") {
				requireAll("pcaPP")
				.PC <- PCAgrid(t(DList[[i]]), k=nPC, method=robust.var, scores=FALSE, center=l1median_HoCr2)
			} else {
				DList[[i]][which(is.na(DList[[i]]))] <- 0 #b/c scale of 0/0
				.PC <- svd(t(DList[[i]]), nu=0) 
				.PC$sdev <- .PC$d/sqrt(ncol(DList[[i]])-1)
				.PC$loadings <- .PC$v
			}
			.totalVar <- c(.totalVar, sum(.PC$sdev^2))
			if(method=="Eigen") {
				.Cov <- cov(t(DList[[i]])) 
				.Cov <- .Cov / (.PC$sdev[1])^2
				.CovSum <- .CovSum + .weight[i] * .Cov
			} else {
				.p <- min(min(sapply(DList, ncol)) - 1, ncol(.PC$loadings))
				.CovSum <- .CovSum + .weight[i] * (.PC$loadings[,1:.p] %*% t(.PC$loadings[,1:.p])) 
			}
		}
		
		if (!all(is.finite(.CovSum))) { #in some cases, Robust method has weird result.
			return(NA)
		}
		
		.eig <- eigen(.CovSum, symmetric=TRUE)
		
		.res <- list()
		.res$x <- lapply(1:length(DList), function(i) list(coord=t(DList[[i]]) %*% .eig$vectors[,1:nPC], totalVar=.totalVar[i]))
		names(.res$x) <- names(DList)
		.res$v <- .eig$vectors[,1:nPC]
		
		if(.scaleAdjust) {
			SDs <- sapply(.res$x, function(x) colSds(x$coord))
			SDsM <- rowMeans(SDs)
			for(i in 1:length(.res$x)) {
				.res$x[[i]]$coord <- sweep(.res$x[[i]]$coord, 2, SDsM/SDs[,i], FUN="*")
			}
		}
	}
	
	return(.res)
}

.MetaSparsePCA <- function(DListF, nPC, maxFeatures=NULL, lambda=NULL, .scale=FALSE, .scaleAdjust=FALSE, max.iter=100, eps=1e-3) {
	union.rec <- function(.list, ...){
		if(length(.list)==1) return(.list[[1]])
		Recall(c(list(union(.list[[1]], .list[[2]], ...)), .list[-(1:2)]), ...)
	}
	features <- union.rec(lapply(DListF,rownames))
	
	stopifnot(maxFeatures <= length(features))
	
	if(is.null(lambda))
		lambda <- rep(length(DListF)/sqrt(length(features)), nPC)
	else
		lambda <- lambda * length(DListF) #to make lambda be comparable to usual spca 
				
	stopifnot(nPC==length(lambda) & nPC <= min(sapply(DListF, ncol)) - 1)
	
	rotation <- matrix(0, length(features), nPC)
	dimnames(rotation) <- list(features, paste("PC",1:ncol(rotation),sep=""))
	
	.DList <- foreach(it=iter(DListF)) %do% { #extended & transposed
		it <- t(scale(t(it), scale=.scale)) 
		.d <- matrix(0, ncol(it), length(features))
		colnames(.d) <- features; rownames(.d) <- colnames(it)
		.d[, match(rownames(it), features)] <- t(it)
		return(.d)
	}
	names(.DList) <- names(DListF)
	..DList <- .DList
	
	#Only consider the first Meta-PC for each iteration
	for(k in 1:nPC) {
		#Beta is the target for Meta-PC 
		#Initialize Betas
		Betas <- foreach(i=1:length(..DList), .combine=cbind) %do% {
			.alpha <- svd(..DList[[i]], nu=0)$v[,1,drop=FALSE]
			.beta <- drop(t(..DList[[i]]) %*% ..DList[[i]] %*% .alpha)
			.norm <- sqrt(sum(.beta^2))
			if(.norm==0) .norm <- 1
			.beta / .norm
		}		
		Betas <- NormalizeMatrix(MetaSoftThreshold(Betas, maxFeatures, lambda[k]))
		
		it <- max.iter
		diff <- 1
		#Find convergence of Betas
		while(it > 1 & diff > eps) {
			BetasNew <- foreach(i=1:length(..DList), .combine=cbind) %do% {
				z <- svd(t(..DList[[i]]) %*% ..DList[[i]] %*% Betas[,i])
				.alpha <- z$u %*% t(z$v)
				
				.beta <- drop(t(..DList[[i]]) %*% ..DList[[i]] %*% .alpha)
				.norm <- sqrt(sum(.beta^2))
				if(.norm==0) .norm <- 1
				.beta / .norm
			}		
			BetasNew <- NormalizeMatrix(MetaSoftThreshold(BetasNew, maxFeatures, lambda[k]))
			
			diff <- max(abs(BetasNew - Betas))
			Betas <- BetasNew
			it <- it - 1
			
			if(it==0 & diff > eps)
				warnings(paste("After", max.iter, "iteration, still not converged!"))
		}
		
		.chosen <- which(rowSums(abs(Betas))>0)
		if(length(.chosen)<1)
			stop("lambda is too large. No features are left")
		#Find Meta-PC which minimize in-between angles
		.CovSum <- Reduce('+', lapply(1:length(..DList), function(i) Betas[.chosen,i,drop=FALSE] %*% t(Betas[.chosen,i,drop=FALSE])))
		.eig <- eigen(.CovSum, symmetric=TRUE)
		
		rotation[.chosen,k] <- .eig$vectors[,1] 
		
		#Orthogonal Projection => however, SparsePCA is usually not orthogonal
		if(k < nPC) {
			for(i in 1:length(..DList)) {
				..DList[[i]] <- ..DList[[i]] %*% (diag(nrow(rotation)) - rotation[,k]%*%t(rotation[,k]))
			}
		}
	}
	
	.res <- list()
	.res$x <- lapply(1:length(.DList), function(i) list(coord=.DList[[i]] %*% rotation[,1:nPC], totalVar=NULL))
	names(.res$x) <- names(.DList)
	.res$v <- rotation
	
	if(.scaleAdjust) {
		SDs <- sapply(.res$x, function(x) colSds(x$coord))
		SDsM <- rowMeans(SDs)
		for(i in 1:length(.res$x)) {
			.res$x[[i]]$coord <- sweep(.res$x[[i]]$coord, 2, SDsM/SDs[,i], FUN="*")
		}
	}
	
	return(.res)
}

