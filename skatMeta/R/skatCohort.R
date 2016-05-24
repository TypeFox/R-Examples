skatCohort <-
function(Z, formula, family = gaussian(), SNPInfo=NULL, snpNames = "Name", aggregateBy = "gene", data=parent.frame(), verbose = FALSE){
	#fit Null model
	if(is.null(SNPInfo)){ 
		warning("No SNP Info file provided: loading the Illumina HumanExome BeadChip. See ?SNPInfo for more details")
		load(paste(find.package("skatMeta"), "data", "SNPInfo.rda",sep = "/"))
		aggregateBy = "SKATgene"
	}
	
	nullmodel <- glm(formula=formula, family = family, data=data)
	res <- residuals(nullmodel, type = "response")
	X1 <- model.matrix(nullmodel)
	n <- nrow(X1)


	env <- environment()
	##check format:
	invisible(check_format_skat(Z, SNPInfo, nullmodel,aggregateBy, snpNames))
	
	##match snps in Z with master list in SNPInfo file 
	mysnps <- colnames(Z)
	
	SNPInfo[,aggregateBy] <- as.character(SNPInfo[,aggregateBy])
	
	SItoZ <- which(colnames(Z) %in% SNPInfo[,snpNames])
	ZtoSI <- which(SNPInfo[,snpNames] %in% colnames(Z))

	which.snps.Z <- colnames(Z) %in% SNPInfo[,snpNames]
	ZtoSI <- match(SNPInfo[,snpNames], mysnps[which.snps.Z])
	#which.snps <- match(mysnps[which.snps.Z],SNPInfo[,snpNames])
	#which.snps <- which(SNPInfo[,snpNames] %in% mysnps[which.snps.Z])
	
	nsnps <- sum(!is.na(ZtoSI)) #sum(which.snps.Z)
	if(nsnps == 0){ 
		stop("no column names in Z match SNP names in the SNP Info file!")
	}
	
	if(verbose){
    	cat("\n Scoring... Progress:\n")
    	pb <- txtProgressBar(min = 0, max = nsnps, style = 3)
    	pb.i <- 0
    }
	
	##fit individual betas/se's
	maf0 <- colMeans(Z,na.rm=TRUE)[which.snps.Z]/2
  	maf0[is.nan(maf0)] <- -1
  
	maf <- maf0[ZtoSI]
	names(maf) <- SNPInfo[,snpNames]
	
	scores <- apply(Z[,which.snps.Z, drop = FALSE],2,function(z){
		if(any(is.na(z))){
      if(all(is.na(z))) z <- rep(0,length(z))
      mz <- mean(z, na.rm=TRUE)
			z[is.na(z)] <- mz
		}
        if (verbose){
				assign("pb.i", get("pb.i",env)+1,env)
				if(get("pb.i", env)%%ceiling(nsnps/100) == 0) setTxtProgressBar(get("pb",env),get("pb.i",env))
			}
		sum(res*z)
		})[ZtoSI]
	scores[is.na(scores)] <- 0
	names(scores) <- SNPInfo[,snpNames]
	
	if(verbose) close(pb)
	
	
	#deal with monomorphic SNPs
	scores[maf == 0] <- 0
	
	#differentiate missing from monomorphic:
	maf[!(SNPInfo[,snpNames] %in% colnames(Z))] <- -1

	#split into genes
	scores 	<- 	split(scores, SNPInfo[,aggregateBy])
	maf 	<- 	split(maf, SNPInfo[,aggregateBy])
	
	##get matrices for projection
	X1 <- sqrt(nullmodel$family$var(nullmodel$fitted))*X1
	AX1 <- solve(crossprod(X1))%*%t(X1)
	
	ngenes <- length(unique(SNPInfo[,aggregateBy]))
	if(verbose){
    	cat("\n Calculating covariance... Progress:\n")
    	pb <- txtProgressBar(min = 0, max = ngenes, style = 3)
    	pb.i <- 0
    }
	##get covariance matrices:
	re <- tapply(SNPInfo[,snpNames], SNPInfo[,aggregateBy],function(snp.names){
		inds <- match(snp.names,colnames(Z))
		mcov <- matrix(0,length(snp.names),length(snp.names))
		if(length(na.omit(inds)) > 0){
			Z0 <- sqrt(nullmodel$family$var(nullmodel$fitted))*as.matrix(Z[,na.omit(inds),drop=FALSE])
			if(any(is.na(Z0))) Z0 <- apply(Z0,2,function(z){
			  if(all(is.na(z))) z <- rep(0,length(z))
       			mz <- mean(z, na.rm=TRUE)
				z[is.na(z)] <- mz
				z
			})
			mcov[!is.na(inds), !is.na(inds)] <- crossprod(Z0) - (t(Z0)%*%X1)%*%(AX1%*%Z0)
		}
		rownames(mcov) <- colnames(mcov) <- snp.names
		if(verbose){
				assign("pb.i", get("pb.i",env)+1,env)
				if(get("pb.i", env)%%ceiling(ngenes/100) == 0) setTxtProgressBar(get("pb",env),get("pb.i",env))		  
		}
		return(mcov)
	},simplify = FALSE)
	sey = sqrt(var(residuals(nullmodel,"pearson"))*(nrow(X1)-1)/(nrow(X1)-ncol(X1)) )
	if(family$family == "binomial") sey = 1
	
	##aggregate
	for(k in 1:length(re)){
		re[[k]] <- list("scores" = scores[[k]], "cov" = re[[k]], "n" =n, "maf" = maf[[k]], "sey" = sey ) 
	}
	if(verbose) close(pb)
	
	attr(re,"family") <-  family$family
	class(re) <- "skatCohort"
	return(re)
}
