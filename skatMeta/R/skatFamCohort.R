skatFamCohort <- function(Z, formula, SNPInfo=NULL, snpNames = "Name", aggregateBy = "gene", data=parent.frame(), fullkins, sparse = TRUE, verbose = FALSE){
	require(coxme)
	##
	if(is.null(SNPInfo)){ 
		warning("No SNP Info file provided: loading the Illumina HumanExome BeadChip. See ?SNPInfo for more details")
		load(paste(find.package("skatMeta"), "data", "SNPInfo.rda",sep = "/"))
		aggregateBy = "SKATgene"
	}
	#if(is.null(id)) id = 1:nrow(data)
	#n <- length(id)
	n <- dim(fullkins)[1]
	
	#tmpidx<-!is.na(match(dimnames(fullkins)[[1]], id))
   # kins<-fullkins[tmpidx, tmpidx]
    #tmpidx<-match(id, dimnames(kins)[[1]])
   # kins<-as.matrix(kins)[tmpidx, tmpidx]
	SNPInfo[,aggregateBy] <- as.character(SNPInfo[,aggregateBy])
	
	kins <- Matrix(kins,sparse=TRUE)
	if(sparse){
		kins[kins < 2 * 2^{-6}] <- 0
		kins <- forceSymmetric(kins)
	}
	
	#fit Null model
	nullmodel <- coxme:::lmekin(formula=update(formula, '~.+ (1|id)'), data=data, varlist = 2*kins)
	nullmodel$theta <- c(nullmodel$vcoef$id*nullmodel$sigma^2,nullmodel$sigma^2)
	
	SIGMA <- nullmodel$theta[1]*2*kins+nullmodel$theta[2]*Diagonal(n)
	X1 <- model.matrix(lm(formula,data=data))
	
	s2 <- sum(nullmodel$theta)
	Om_i <- solve(SIGMA/s2)
	
	#rotate data:
	res <- as.vector(nullmodel$res)	
	
	##match snps in Z with master list in SNPInfo file 
	mysnps <- colnames(Z)
	
	which.snps.Z <- colnames(Z) %in% SNPInfo[,snpNames]
	#which.snps <- match(mysnps[which.snps.Z],SNPInfo[,snpNames])
	ZtoSI <- match(SNPInfo[,snpNames], mysnps[which.snps.Z])

	##fit individual betas/se's
	maf0 <- colMeans(Z,na.rm=TRUE)[which.snps.Z]/2
  	maf0[is.nan(maf0)] <- -1

	maf <- maf0[ZtoSI]
	names(maf) <- SNPInfo[,snpNames]

	nsnps <- sum(!is.na(ZtoSI))
	
	if(nsnps == 0){ 
		stop("no column names in Z match SNP names in the SNP Info file!")
	}
	
	env <- environment()
	if(verbose){
    	cat("\n Scoring... Progress:\n")
    	pb <- txtProgressBar(min = 0, max = nsnps, style = 3)
    	pb.i <- 0
    }

	trOmi <- t(res) * s2 / nullmodel$theta[2]
	scores <- apply(Z[,which.snps.Z,drop=FALSE],2,function(z){
		if(any(is.na(z))){
		  if(all(is.na(z))) z <- rep(0,length(z))
			mz <- mean(z, na.rm=TRUE)
			z[is.na(z)] <- mz
		}
		if (verbose){
				assign("pb.i", get("pb.i",env)+1,env)
				if(get("pb.i", env)%%ceiling(nsnps/100) == 0) setTxtProgressBar(get("pb",env),get("pb.i",env))
		}
		as.numeric(trOmi%*%z)
		})[ZtoSI]
	
	scores[is.na(scores)] <- 0
	names(scores) <- SNPInfo[,snpNames]

	#deal with monomorphic SNPs
	scores[maf == 0] <- 0

	#differentiate missing from monomorphic:
	maf[!(SNPInfo[,snpNames] %in% colnames(Z))] <- -1

	#split into genes
	scores <- split(scores, SNPInfo[,aggregateBy])
	maf <- split(maf, SNPInfo[,aggregateBy])
	
	##get matrices for projection
	AX1 <- solve(t(X1)%*%Om_i%*%X1)%*%t(X1)%*%Om_i
	
	##get covariance matrices:
	if(verbose) close(pb)

	ngenes <- length(unique(SNPInfo[,aggregateBy]))
	if(verbose){
    	cat("\n Calculating covariance... Progress:\n")
    	pb <- txtProgressBar(min = 0, max = ngenes, style = 3)
    	pb.i <- 0
    }	
	re <- tapply(SNPInfo[,snpNames], SNPInfo[,aggregateBy],function(snp.names){
		inds <- match(snp.names,colnames(Z))
		mcov <- matrix(0,length(snp.names),length(snp.names))
		rownames(mcov) <- colnames(mcov) <- snp.names
		if(length(na.omit(inds)) > 0){
			Z0 <- as.matrix(Z[,na.omit(inds), drop = FALSE])
			if(any(is.na(Z0))) Z0 <- apply(Z0,2,function(z){
			  if(all(is.na(z))) z <- rep(0,length(z))
        mz <- mean(z, na.rm=TRUE)
				z[is.na(z)] <- mz
				z
			})
			mcov[!is.na(inds), !is.na(inds)] <- as.matrix(t(Z0)%*%Om_i%*%Z0 - (t(Z0)%*%Om_i%*%X1)%*%(AX1%*%Z0))
		}
		if(verbose){
			assign("pb.i", get("pb.i",env)+1,env)
			if(get("pb.i", env)%%ceiling(ngenes/100) == 0) setTxtProgressBar(get("pb",env),get("pb.i",env))		  
		}

		return(mcov)
	},simplify=FALSE)
	if(verbose) close(pb)

	##aggregate
	sey = sqrt(s2)
	for(k in 1:length(re)){
		re[[k]] <- list("scores" = scores[[k]],"cov" = as.matrix(re[[k]]), "n" =n, "maf" = maf[[k]], "sey" = sey) 
	}

	class(re) <- "skatCohort"
	return(re)
}