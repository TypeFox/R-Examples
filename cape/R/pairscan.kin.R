pairscan.kin <-
function(data.obj, geno.obj = NULL, covar, scan.what, marker.pairs, kin.full.geno, sample.kinship, num.kin.samples, n.per.sample, verbose = TRUE, n.cores = 2){

	m = NULL #for appeasing R CMD check
	

	pheno <- get.pheno(data.obj, scan.what)	
	num.pheno <- dim(pheno)[2]
	results.list <- vector(mode = "list", length = num.pheno)
	names(results.list) <- colnames(pheno)

	covar.info <- get.covar(data.obj, covar)
	is.char <- as.logical(is.na(suppressWarnings(as.numeric(marker.pairs[1,1]))))
	if(is.char){
		covar.names <- get.marker.name(data.obj, covar.info$covar.names)
		}else{
		covar.names <- get.marker.num(data.obj, covar.info$covar.names)	
		}
	num.covar <- length(covar.names)
	p.covar.locale <- which(covar.info$covar.type == "p")
	num.p.covar <- length(p.covar.locale)
	#============================================================================
	#internal functions
	#============================================================================

	get.marker.pair.stats <- function(m, new.geno, new.pheno, new.covar, err.cov){
		int.term = matrix(solve(err.cov) %*% new.geno[,m[1]]*new.geno[,m[2]], ncol = 1)
		pairscan.results <- one.pairscan.parallel(data.obj, phenotype.vector = new.pheno, genotype.matrix = new.geno, int = int.term, covar.vector = new.covar, paired.markers = matrix(m, ncol = 2), n.perm = 0, verbose = FALSE, run.parallel = FALSE, n.cores = n.cores)
		if(is.null(pairscan.results[[1]])){
			marker.num <- get.marker.num(data.obj, m)
			dummyV <- c(marker.num, rep(NA, 3))
			results <- list("effects" = dummyV, "se" = dummyV, "cov" = c(dummyV, rep(NA, 4)))
			}else{
			results <- list("effects" = pairscan.results[[1]]$pairscan.effects, "se" = pairscan.results[[1]]$pairscan.se, "cov" = pairscan.results[[1]]$model.covariance)
			}
		return(results)
		}


	get.covar.stats <- function(m, new.geno, new.pheno, new.covar, err.cov){
		covar.locale <- matrix(which(covar.names == m[2]), ncol = 1)
		int.term = solve(err.cov) %*% new.geno[,m[1]]*new.covar[,covar.locale]
		#the genotype.matrix does not have the phenotypic covariates in
		#it, so append these
		pairscan.results <- one.pairscan.parallel(data.obj, phenotype.vector = new.pheno, genotype.matrix = cbind(new.geno, new.covar[,covar.locale,drop=FALSE]), int = int.term, covar.vector = new.covar[,-covar.locale], paired.markers = matrix(m, ncol = 2), n.perm = 0, verbose = FALSE, run.parallel = FALSE, n.cores = n.cores)
		if(is.null(pairscan.results[[1]])){
			marker.num <- get.marker.num(data.obj, m)
			dummyV <- c(marker.num, rep(NA, 3))
			results <- list("effects" = dummyV, "se" = dummyV, "cov" = c(dummyV, rep(NA, 4)))
			}else{
			results <- list("effects" = pairscan.results[[1]]$pairscan.effects, "se" = pairscan.results[[1]]$pairscan.se, "cov" = pairscan.results[[1]]$model.covariance)
			}
		return(results)
		}
	#============================================================================
	

	#calculate kinship matrices
	#The genotype matrix stored in the geno.obj can be used
	#preferrentially, even if only a subset of markers are
	#being scanned. In this case, take the genotype matrix
	#directly from the genotype object and filter it for the
	#correct individuals. If there is no genotype object, we
	#use the genotype matrix in the data object. If the full
	#genotype matrix is not requested, use get.geno() to use
	#the markers that are being scanned in the kinship matrices
	if(kin.full.geno && !is.null(geno.obj)){
		full.geno <- geno.obj$geno
		#match up the individuals
		ind.locale <- match(rownames(pheno), rownames(full.geno))
		full.geno <- full.geno[ind.locale,]
		}else{
		#do not use the phenotypic covariates in calculating the kinship matrices
		full.geno <- get.geno.with.covar(data.obj, geno.obj, g.covar = TRUE, p.covar = FALSE, for.pairscan = FALSE)
		}
		
		used.markers <- colnames(full.geno)
		if(kin.full.geno){		
			used.chr <- get.marker.chr(geno.obj, used.markers)			
			}else{
			used.chr <- get.marker.chr(data.obj, used.markers)				
			}

		u_used.chr <- unique(used.chr)
		testing.covar <- as.logical(num.p.covar > 0)
		
		kin.mats <- kinship(geno = full.geno, all.chr = used.chr, chr.which = u_used.chr, chr.pairs = TRUE, sample.kinship = sample.kinship, num.samples = num.kin.samples, n.per.sample = n.per.sample, add.full.kin = testing.covar, verbose = verbose, n.cores = n.cores)


	#get all the chromosome pairs. We only want to calculate the 
	#kin-corrected phenotype and genotype once for each chromosome
	#pair, so we calculate stats per chromosome pair
	chr.pairs <- pair.matrix(u_used.chr, self.pairs = TRUE)		
	#determine which chromosome each marker is on
	
	if(kin.full.geno || is.null(geno.obj)){
		marker.chr <- apply(marker.pairs, 2, function(x) get.marker.chr(data.obj, x))		
		}else{
		marker.chr <- apply(marker.pairs, 2, function(x) get.marker.chr(geno.obj, x))
		}
	
	
	for(p in 1:num.pheno){
		if(verbose){
			cat("\nScanning phenotype ", colnames(pheno)[p], ":\n", sep = "")
			}

		covar.vector <- covar.info$covar.table
		pheno.vector <- pheno[,p]

		start.ind <- 1

		effects.mat <- matrix(NA, nrow = dim(marker.pairs)[1], ncol = 5)
		se.mat <- matrix(NA, nrow = dim(marker.pairs)[1], ncol = 5)
		cov.mat <- matrix(NA, nrow = dim(marker.pairs)[1], ncol = 9)

		for(chp in 1:dim(chr.pairs)[1]){
		# for(chp in 1:25){	
			if(verbose){cat("\tChromosome(s)", unique(chr.pairs[chp,]), "...\n")}

			#find all the marker pairs between the chromosomes
			chr.markers <- get.chr.markers(chr.pair = chr.pairs[chp,], marker.chr, paired.markers = marker.pairs)
			num.chr.pairs <- dim(chr.markers)[1]
			
			#If there are markers on this chromosome pair
			if(num.chr.pairs > 0){
				#calculate the corrected genotype and phenotype for this chromosome pair
				kin.dat <- kinship.on.the.fly(kin.mats, geno = full.geno, chr1 = chr.pairs[chp,1], chr2 = chr.pairs[chp,2], phenoV = pheno.vector, covarV = covar.vector)


				new.pheno <- kin.dat$corrected.pheno
				new.geno <- kin.dat$corrected.geno
				new.covar <- kin.dat$corrected.covar
				err.cov <- kin.dat$err.cov
			
				cl <- makeCluster(n.cores)
				registerDoParallel(cl)
				pairscan.results <- foreach(m = t(chr.markers)) %dopar% {
					get.marker.pair.stats(m, new.geno, new.pheno, new.covar, err.cov)
					}
				stopCluster(cl)
				effects.mat[start.ind:(start.ind+num.chr.pairs-1),] <- matrix(unlist(lapply(pairscan.results, function(x) x$effects)), nrow = num.chr.pairs, byrow = TRUE)
				se.mat[start.ind:(start.ind+num.chr.pairs-1),] <- matrix(unlist(lapply(pairscan.results, function(x) x$se)), nrow = num.chr.pairs, byrow = TRUE)
				cov.mat[start.ind:(start.ind+num.chr.pairs-1),] <- matrix(unlist(lapply(pairscan.results, function(x) x$cov)), nrow = num.chr.pairs, byrow = TRUE)
			
				start.ind <- start.ind+num.chr.pairs
				} #end case for there being markers on this chromosome pair
			} #end looping through chromosome pairs
			
			#if there are covariates to test explicitly	
			if(num.p.covar > 0){
				for(cv in 1:num.p.covar){
				if(verbose){cat("\tCovariate:", covar.names[p.covar.locale[cv]], "\n")}
				cv.marker.locale <- c(which(marker.pairs[,1] == covar.names[p.covar.locale[cv]]), which(marker.pairs[,2] == covar.names[p.covar.locale[cv]]))
				cv.markers <- marker.pairs[cv.marker.locale,,drop=FALSE]
				num.cv.pairs <- dim(cv.markers)[1]
				
				if(num.cv.pairs > 0){
					#calculate the corrected genotype and phenotype for this chromosome pair
					kin.dat <- kinship.on.the.fly(kin.mats, geno = full.geno, chr1 = NULL, chr2 = NULL, phenoV = pheno.vector, covarV = covar.vector)

					new.pheno <- kin.dat$corrected.pheno
					new.geno <- kin.dat$corrected.geno
					new.covar <- kin.dat$corrected.covar
					if(is.char){
						colnames(new.covar) <- get.marker.name(data.obj, as.numeric(colnames(new.covar)))
						}
					err.cov <- kin.dat$err.cov
				
					cl <- makeCluster(n.cores)
					registerDoParallel(cl)
					pairscan.results <- foreach(m = t(cv.markers)) %dopar% {
						get.covar.stats(m, new.geno, new.pheno, new.covar, err.cov)
						}
					stopCluster(cl)
					
					effects.mat[start.ind:(start.ind+num.cv.pairs-1),] <- matrix(unlist(lapply(pairscan.results, function(x) x$effects)), nrow = num.cv.pairs, byrow = TRUE)
					se.mat[start.ind:(start.ind+num.cv.pairs-1),] <- matrix(unlist(lapply(pairscan.results, function(x) x$se)), nrow = num.cv.pairs, byrow = TRUE)
					cov.mat[start.ind:(start.ind+num.cv.pairs-1),] <- matrix(unlist(lapply(pairscan.results, function(x) x$cov)), nrow = num.cv.pairs, byrow = TRUE)
				
					start.ind <- start.ind+num.chr.pairs
					}#end case for when there are pairs of markers with covariates
				} #end looping through markers paired with covariates
			}#end case for when there are phenotypic covariates to test
		
			colnames(effects.mat) <- colnames(pairscan.results[[1]]$effects)
			colnames(se.mat) <- colnames(pairscan.results[[1]]$se)
				
			effects.mat <- apply(effects.mat, 2, as.numeric)
			se.mat <- apply(se.mat, 2, as.numeric)
									
			#add the results for the phenotype
			pheno.results <- list(effects.mat, se.mat, cov.mat)
			names(pheno.results) <- c("pairscan.effects", "pairscan.se", "model.covariance")
			results.list[[p]] <- pheno.results
		}	#end looping over phenotypes	
 

		return(results.list)


}
