singlescan <-
function(data.obj, geno.obj = NULL, n.perm = NULL, covar = NULL, alpha = c(0.01, 0.05), scan.what = c("eigentraits", "raw.traits"), use.kinship = FALSE, kin.full.geno = FALSE, run.parallel = TRUE, sample.kinship = FALSE, num.kin.samples = 1000, n.per.sample = 10, verbose = FALSE, overwrite.alert = TRUE, n.cores = 2) {

	if(overwrite.alert){
		choice <- readline(prompt = "Please make sure you assign the output of this function to a singlescan.obj, and NOT the data.obj. It will overwrite the data.obj.\nDo you want to continue (y/n) ")
		if(choice == "n"){stop()}
		}


	gene <- get.geno(data.obj, geno.obj)
	pheno <- get.pheno(data.obj, scan.what)			
	n.phe = dim(pheno)[2]
	chr.which <- unique(data.obj$chromosome)
	
	#get the covariates and assign the variables
	#to the local environment
	covar.info <- get.covar(data.obj, covar)
	for(i in 1:length(covar.info)){
		assign(names(covar.info)[i], covar.info[[i]])
		}
		
		
	if(is.null(covar) && !is.null(covar.names)){
		choice <- readline(prompt = "There are covariates in the data object, but none are specified. Would you like to continue? (y/n)")
		if(choice == "n"){stop()}
		}
		
	if(is.null(covar)){
		covar.names = "none"
		covar.loc = NULL
		covar.table = NULL
		}

	#if we are using covariates, we need to calculate
	#and extra kinship matrix using the full genome
	if(!is.null(covar)){
		add.full.kin <- TRUE
		}else{
		add.full.kin <- FALSE
		}
	
	if(is.null(n.perm) || n.perm < 2){alpha = "none"}
	singlescan.obj <- vector(mode = "list", length = 4)
	names(singlescan.obj) <- c("alpha", "alpha.thresh", "covar", "singlescan.results")
	singlescan.obj$covar <- colnames(covar.table)
	singlescan.obj$alpha <- alpha

	#==================================================================
	#if we are using a kinship correction, make sure the phenotypes
	#are mean-centered and there are no missing values in the genotype
	#matrix.
	#==================================================================
	if(use.kinship){
		pheno.means <- apply(pheno, 2, mean)
		tol = 0.01
		non.zero <- intersect(which(pheno.means > 0+tol), which(pheno.means < 0-tol))
		if(length(non.zero) > 0){
			warning("Phenotypes must be mean-centered before performing kinship corrections.")
			cat("Mean-centering phenotypes using norm.pheno()")
			data.obj <- norm.pheno(data.obj)
			# if(length(type.choice) > 0){ #if we find any, use the eigentrait matrix
				# pheno <- data.obj$ET
				# }else{
				# pheno <- data.obj$pheno #otherwise, use the raw phenotype matrix
				# }
			}
		missing.vals <- which(is.na(gene))
			if(length(missing.vals) > 0){
				stop("There are missing values in the genotype matrix. Please use impute.geno().")
				}
			}
	#==================================================================


	#==================================================================
	#first do the permutations to get the significance threshold
	#results will be compared to the significance threshold to 
	#determine	
	#==================================================================
	if(!is.null(n.perm) && n.perm >= 2){
		if(verbose){
		cat("\nPerforming permutations to calculate significance threshold...\n")
		}
		if(run.parallel){
			alpha.thresh <- genome.wide.threshold.1D.parallel(data.obj, geno.mat = gene, n.perm = n.perm, alpha = alpha, scan.what = scan.what, verbose = verbose, n.cores = n.cores)
			}else{
			alpha.thresh <- genome.wide.threshold.1D(data.obj, geno.mat = gene, n.perm = n.perm, alpha = alpha, scan.what = scan.what, verbose = verbose)
			}
	singlescan.obj$alpha.thresh <- alpha.thresh
		}else{
		if(verbose){message("Permutations are not being calculated\n\tTo calculate permutations n.perm must be greater than 2.\n")}	
		singlescan.obj$alpha.thresh = "no permutations"
		}
	#==================================================================
	

		#===========================================================================
		#internal functions
		#===========================================================================		

		#This function gets regression statistics with a
		#covariate table
		get.stats <- function(phenoV, markerV, covarV = NULL){
			
			marker.var <- var(markerV, na.rm = TRUE)
			if(marker.var == 0){
				return(rep(NA, 4))
				}

			if(use.kinship){
				if(!is.null(covarV)){
					model <- lm(phenoV~0+cbind(covarV, markerV))#remove intercept from model	
					}else{
					model <- lm(phenoV~0+markerV)#remove intercept from model		
					}
				}else{
				if(!is.null(covarV)){
					model <- lm(phenoV~cbind(covarV, markerV))
					}else{
					model <- lm(phenoV~markerV)
					}
				}
			#take the last line of coefficients.
			model.coef <- summary(model)$coefficients
			slope <- model.coef[dim(model.coef)[1],1]
			se <- model.coef[dim(model.coef)[1],2]
			t.stat <- abs(model.coef[dim(model.coef)[1],3])
			p.val <- model.coef[dim(model.coef)[1],4]
							
			#put together all the statistics we want to keep
			#we keep the absolute value of the t statistic,
			#the p value, and the covariate flag
			
			table.row <- c(slope, se, t.stat, p.val)
			rm("model.coef", "slope", "se", "t.stat", "p.val", "model")
			return(table.row)
			}		
		
		get.chr.stats <- function(ch, phenoV, covar.mat){
			if(verbose){cat("\tCh", ch, "...", sep = "")}
			marker.locale <- which(data.obj$chromosome == ch)
			if(use.kinship){
				kin.obj <- kinship.on.the.fly(kin.mats = kin.mats, geno = kin.gene.calc, chr1 = ch, chr2 = ch, phenoV = phenoV, covarV = covar.mat)
				cor.geno <- kin.obj$corrected.geno
				cor.pheno <- kin.obj$corrected.pheno
				cor.covar <- kin.obj$corrected.covar
				}else{
				cor.geno <- gene
				cor.pheno <- phenoV
				cor.covar <- covar.mat
				}
				
				#calculate the singlescan for the chr markers using the kinship object
				if(run.parallel){
					just.markers <- cor.geno[,marker.locale,drop=FALSE]
					cl <- makeCluster(n.cores)
					registerDoParallel(cl)
					marker.stats <- foreach(m = just.markers, .combine = "cbind", .export = c("get.stats")) %dopar% {
						get.stats(phenoV = cor.pheno, markerV = m, covarV = cor.covar)
						}
					stopCluster(cl)
					}else{
	
				marker.stats <- matrix(NA, ncol = length(marker.locale), nrow = 4)
				for(m in 1:length(marker.locale)){
					report.progress(m, length(marker.locale))
					marker.stats[,m] <- get.stats(phenoV = cor.pheno, markerV = cor.geno[,marker.locale[m]], covarV = cor.covar)
					}
				if(verbose){cat("\n")}
		
					}
				if(length(marker.locale) == 1){
					marker.stats <- matrix(marker.stats, ncol = 1)
					}
				colnames(marker.stats) <- data.obj$marker.num[marker.locale]
				return(marker.stats)
				}
		
		
		#===========================================================================
		
		pheno.stats <- vector(mode = "list", length = n.phe)
	
		if(use.kinship){
	
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
				kin.gene.calc <- geno.obj$geno
				#match up the individuals
				ind.locale <- match(rownames(pheno), rownames(full.geno))
				full.geno <- full.geno[ind.locale,]
				geno.chr.calc <- get.marker.chr(data.obj, colnames(kin.gene.calc))
				}else{
				kin.gene.calc <- get.geno.with.covar(data.obj, geno.obj, p.covar = FALSE, for.pairscan = FALSE)
				geno.chr.calc <- get.marker.chr(data.obj, colnames(kin.gene.calc))
				}

			chr.which <- unique(geno.chr.calc)
			
			#get the kinship matrices for all LOO chromosomes
			kin.mats <- kinship(geno = kin.gene.calc, all.chr = geno.chr.calc, chr.pairs = FALSE, sample.kinship = sample.kinship, num.samples = num.kin.samples, n.per.sample = n.per.sample, run.parallel = run.parallel, add.full.kin = add.full.kin, verbose = verbose, n.cores = n.cores)	
			
			}
		
		for(ph in 1:n.phe){
			if(verbose){cat("\nScanning", colnames(pheno)[ph], "\n")}
			phenoV <- pheno[,ph]
			# if(use.svm){
				# covar.mat <- cbind(covar.table, svm.covar[,ph])
				# }else{
				covar.mat <- covar.table
				# }
			pheno.stats[[ph]] <- lapply(chr.which, function(x) get.chr.stats(x, phenoV, covar.mat))
			}
	
			
		#reorganize and label pheno.stats
		org.pheno.stats <- vector(mode = "list", length = n.phe)
		names(org.pheno.stats) <- colnames(pheno)
		for(ph in 1:length(pheno.stats)){
			num.markers <- sum(unlist(lapply(pheno.stats[[ph]], function(x) dim(x)[2])))
			stats.mat <- matrix(NA, ncol = 4, nrow = num.markers)
			colnames(stats.mat) <- c("slope", "se", "t.stat", "p.val")
			rownames(stats.mat) <- unlist(lapply(pheno.stats[[ph]], function(x) dimnames(x)[[2]]))
			start.pos <- 1
			for(m in 1:length(pheno.stats[[ph]])){
				nrow.chunk <- dim(pheno.stats[[ph]][[m]])[2]
				stats.mat[start.pos:(start.pos+nrow.chunk-1),] <- t(pheno.stats[[ph]][[m]])
				start.pos <- start.pos+nrow.chunk
				}
			org.pheno.stats[[ph]] <- stats.mat
			}

		pheno.stats <- org.pheno.stats
	

		#calculate statistics for the covariates
		if(!is.null(covar)){
			if(verbose){cat("Calculating statistics for covariates.\n")}
			covar.stats <- vector(mode = "list", length = n.phe)
			for(ph in 1:n.phe){ 
				if(use.kinship){
					#calculate the kinship object using the entire genome
					kin.obj <- kinship.on.the.fly(kin.mats = kin.mats, geno = kin.gene.calc, chr1 = NULL, chr2 = NULL, phenoV = pheno[,ph], covarV = covar.table)
					cor.geno <- kin.obj$corrected.geno
					cor.pheno <- kin.obj$corrected.pheno
					cor.covar <- kin.obj$corrected.covar
					}else{
					cor.geno <- gene
					cor.pheno <- pheno[,ph]
					cor.covar <- covar.table
					}
					#calculate the singlescan using the covariates
					covar.stats[[ph]] <- apply(matrix(1:dim(cor.covar)[2], nrow = 1), 2, function(x) get.stats(cor.pheno, cor.covar[,x], NULL))
					}

			for(ph in 1:length(pheno.stats)){
				pheno.stats[[ph]] <- rbind(pheno.stats[[ph]], matrix(covar.stats[[ph]], nrow = dim(covar.table)[2], byrow = TRUE))
				rownames(pheno.stats[[ph]])[which(rownames(pheno.stats[[ph]]) == "")] <- colnames(covar.table)
				}
			} #end case for when we have covariates


					
	
		singlescan.obj$singlescan.results <- pheno.stats
		
		return(singlescan.obj)
	
	}
