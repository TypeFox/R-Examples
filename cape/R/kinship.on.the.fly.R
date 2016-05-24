kinship.on.the.fly <-
function(kin.mats, geno, chr1 = NULL, chr2 = NULL, phenoV, covarV = NULL){

	get.g=function(pair = NULL, phenotype, covarV){
		
		if(is.null(pair)){
			pair.name <- "full.geno"
			}else{
			pair.name <- paste(pair, collapse = ",")
			}
		kin.mat.locale <- which(names(kin.mats) == pair.name)
		K <- kin.mats[[kin.mat.locale]]
		
		#if we are correcting the covariate only don't put it in the model
		if(is.null(covarV) || is.null(pair)){
			model = regress(as.vector(phenotype)~1,~K, pos = c(TRUE, TRUE))	
			}else{
			model = regress(as.vector(phenotype)~covarV, ~K, pos = c(TRUE, TRUE))		
			}
				
		#This err.cov is the same as err.cov in Dan's code using estVC
		err.cov = summary(model)$sigma[1]*K+summary(model)$sigma[2]*diag(nrow(K))

		
		eW = eigen(err.cov, symmetric = TRUE)
	    if(min(eW$values) < 0 && abs(min(eW$values)) > sqrt(.Machine$double.eps)){
		      }else{
	      	eW$values[eW$values <= 0] = Inf
	      	} 
	      	err.cov = eW$vector %*% diag(eW$values^-0.5) %*% t(eW$vector)

		new.pheno <- err.cov %*% phenotype
		new.geno <- err.cov %*% geno; #hist(new.geno)
		rownames(new.geno) <- rownames(geno)
		if(!is.null(covarV)){
			new.covar <- err.cov %*% covarV
			}else{
			new.covar <- NULL	
			}
		
		results = list(err.cov, new.pheno, new.geno, new.covar)
		names(results) <- c("err.cov", "corrected.pheno", "corrected.geno", "corrected.covar")
		return(results)
		}

	result <- get.g(pair = c(chr1, chr2), phenotype = phenoV, covarV = covarV)
	
	return(result)
	
}
