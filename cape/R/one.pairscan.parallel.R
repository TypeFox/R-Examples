one.pairscan.parallel <-
function(data.obj, phenotype.vector, genotype.matrix, int = NULL, covar.vector = NULL, paired.markers, n.perm = 0, run.parallel = TRUE, verbose = FALSE, n.cores = 2){
		
		m = p = NULL #for appeasing R CMD check

		is.char <- as.logical(is.na(suppressWarnings(as.numeric(paired.markers[1,1]))))
		num.markers <- dim(genotype.matrix)[2]
		if(length(covar.vector) > 0){
		if(is.char){
			covar.names <- get.marker.name(data.obj, colnames(covar.vector))
			}else{
			covar.names <- get.marker.num(data.obj, colnames(covar.vector))
			}
			}else{
			covar.names <- NULL	
			}
		#============================================================================
		# check to see that the covariates are not redundant and are linearly independent
		#============================================================================
		use.covars <- as.logical(length(covar.vector) > 0)
		if(use.covars > 0){
			cov.mat <- covar.vector
			
			#remove the NAs and check the matrix for rank
			not.na.locale <- which(!is.na(rowSums(cov.mat)))
			no.na.cov <- cov.mat[not.na.locale,,drop=FALSE]
			
			design.cov <- cbind(rep(1, dim(no.na.cov)[1]), no.na.cov)
			rank.cov <- rankMatrix(design.cov)
			if(rank.cov[[1]] < dim(design.cov)[2]){
				stop("The covariate matrix does not appear to be linearly independent.\nIf you are using dummy variables for groups, leave one of the groups out.")
				}
			
			cor.mat <- cor(no.na.cov)
			diag(cor.mat) <- 0
			perfect.cor <- which(abs(signif(cor.mat, 2)) == 1)
			if(length(perfect.cor) > 0){
				stop("Check the covariates. There appears to be at least one pair of redundant covariates.")
				}
			num.covar <- dim(covar.vector)[2]	
			}else{
			num.covar = 0
			}			
		#============================================================================
		
			
		#============================================================================
		#internal functions
		#============================================================================
		
			
		get.model.results <- function(marker.names, m1, m2, int.term = NULL, testing.covar = FALSE){
				
				#if we are testing a covariate, pull it out of the 
				#covariate.table
				if(testing.covar){
					covar.locale <- which(covar.names %in% marker.names)
					if(length(covar.locale) > 0){
						new.covar.vector <- covar.vector[,-covar.locale,drop=FALSE]
						}else{
						new.covar.vector <- covar.vector	
						}
					}else{
					new.covar.vector <- covar.vector
					}
					
				if(is.null(new.covar.vector) || dim(new.covar.vector)[2] == 0){
					new.covar.vector <- NULL
					}
					
				design.mat <- cbind(rep(1, length(m1)), new.covar.vector, m1, m2)
				
				missing.rows <- which(is.na(rowSums(design.mat)))
				if(length(missing.rows) > 0){
					design.mat <- design.mat[-missing.rows,]
					}
				rank.check <- rankMatrix(design.mat)[[1]]

				if(rank.check < dim(design.mat)[2]){
					return(NULL)
					}

				#Do the linear regression with the covariates, the two markers 
				#individually, and the interaction between the two markers
				#put the covars first so the marker effects come last
				if(!is.null(new.covar.vector)){
					if(is.null(int.term)){
						model <- lm(phenotype.vector ~ new.covar.vector + m1 + m2 + m1:m2, model = FALSE, x = FALSE, y = FALSE, qr = TRUE)
						}else{
						model <- lm(phenotype.vector ~ 0 + new.covar.vector + m1 + m2 + int.term, model = FALSE, x = FALSE, y = FALSE, qr = TRUE)	
						}
					}else{
						if(is.null(int.term)){
							model <- lm(phenotype.vector ~ m1 + m2 + m1:m2, model = FALSE, x = FALSE, y = FALSE, qr = TRUE)	
							}else{
							model <- lm(phenotype.vector ~ 0 + m1 + m2 + int.term, model = FALSE, x = FALSE, y = FALSE, qr = TRUE)		
							}
						}
				model.summ <- summary(model)
				
				#occasionally permutations result in non linearly dependent matrices.
				#if this is the case, return NULL. This triggers the permutation
				#script to generate another permutation.
				if(length(which(is.na(coefficients(model)))) > 0){ 
					return(NULL)
					}
					
				
				#take the last 3 terms
				keep.results <- tail(1:length(coef(model)), 3)
				model.effects <- model.summ$coefficients[keep.results,"Estimate"]
				model.se <- model.summ$coefficients[keep.results,"Std. Error"]

				#calculate the covariance matrix of the model parameters
				model.cov.mat <- vcov(model)
				dim.mat <- dim(model.cov.mat)[1] #get the dimensions of the matrix. We want the last three rows and last three columns (the covariance matrix for m1, m2, and m1:m2)
				cov.mat <- model.cov.mat[(dim.mat-2):dim.mat, (dim.mat-2):dim.mat]
				model.cov.results <- as.vector(cov.mat)

				results <- list(model.effects, model.se, model.cov.results)
				names(results) <- c("model.effects", "model.se", "model.cov")
				return(results)
				}
				
			get.pair.results <- function(m.pair){
				#check the markers for a covariate
				covar.check <- which(m.pair %in% covar.names)
				if(length(covar.check) > 0){
					testing.covar = TRUE
					}else{
					testing.covar = FALSE	
					}

				#get the marker identities
				marker1 <- genotype.matrix[,as.character(m.pair[1])]
				marker2 <- genotype.matrix[,as.character(m.pair[2])]

				if(is.null(int)){
					marker.pair.results <- get.model.results(marker.names = m.pair, m1 = marker1, m2 = marker2, testing.covar = testing.covar)
					}else{
					marker.pair.results <- get.model.results(marker.names = m.pair, m1 = marker1, m2 = marker2, int.term = int, testing.covar = testing.covar)
					}

				return(marker.pair.results)				
				}
				
		one.perm <- function(perm.num){
			m.pair <- paired.markers[sample(1:dim(paired.markers)[1], 1),]
			rnd.order <- sample(1:dim(genotype.matrix)[1], dim(genotype.matrix)[1])
			marker1 <- genotype.matrix[rnd.order,as.character(m.pair[1])]
			marker2 <- genotype.matrix[rnd.order,as.character(m.pair[2])]

			if(is.null(int)){
				marker.pair.results <- get.model.results(marker.names = m.pair, m1 = marker1, m2 = marker2)
				}else{
				marker.pair.results <- get.model.results(marker.names = m.pair, m1 = marker1, m2 = marker2, int.term = int)	
				}
			
			marker.pair.results$"pair.used" <- m.pair
			return(marker.pair.results)				
			}
		#============================================================================
		#end of internal functions
		#============================================================================
		
		# each.iter.time <- rep(NA, n.pairs)
		#we can calculate the number of genes
		#from the input data
		n.pairs <- dim(paired.markers)[1]

		# if(verbose){cat("\tCalculating pair statistics...\n")}
		if(run.parallel){
			cl <- makeCluster(n.cores)
			registerDoParallel(cl)
			pair.results <- foreach(m = t(paired.markers), .export = "phenotype.vector") %dopar% {
				get.pair.results(m)
				}
			stopCluster(cl)
			}else{
			pair.results <- lapply(1:dim(paired.markers)[1], function(x) get.pair.results(paired.markers[x,]))
			}
			
		#Filter out the results with null values
		good.results <- which(unlist(lapply(pair.results, function(x) length(x$model.effects))) > 0)
		if(length(good.results) > 0){
			pair.results <- pair.results[good.results]
			all.model.effects <- matrix(unlist(lapply(pair.results, function(x) x$model.effects)), nrow = length(pair.results), byrow = TRUE)
			all.model.se <- matrix(unlist(lapply(pair.results, function(x) x$model.se)), nrow = length(pair.results), byrow = TRUE)
			all.model.cov <- matrix(unlist(lapply(pair.results, function(x) x$model.cov)), nrow = length(pair.results), byrow = TRUE)
	
		    #assign column names to the results tables
			#the column names represent the possible beta
			#coefficients we can get. The intercept, all
			#possible covariates, marker1 and marker2,
			#and the interaction marker1:marker2
			column.names <- c("marker1", "marker2", "marker1:marker2")	
			colnames(all.model.effects) <- colnames(all.model.se) <- column.names
			
			#add the marker pair names to the results tables and name the columns
	      	paired.marker.num <- matrix(apply(paired.markers[good.results,,drop=FALSE], 2, function(x) get.marker.num(data.obj, x)), ncol = 2, byrow = FALSE)
	      	colnames(paired.marker.num) <- c("marker.name1", "marker.name2")
			final.effects.table <- cbind(paired.marker.num, all.model.effects)
			final.se.table <- cbind(paired.marker.num, all.model.se)
			final.cov.table <- all.model.cov
	
		
			phenotype.results <- list(final.effects.table, final.se.table, final.cov.table)
			names(phenotype.results) <- c("pairscan.effects", "pairscan.se", "model.covariance")
			}else{
			phenotype.results <- NULL	
			}

		
		if(n.perm > 0){
			if(verbose){cat("\tCalculating permutations...\n")}

			if(run.parallel){
				cl <- makeCluster(n.cores)
				registerDoParallel(cl)
				perm.results <- foreach(p = 1:n.perm, .export = "phenotype.vector") %dopar% {
					one.perm(p)
					}
				stopCluster(cl)
				}else{
				perm.results <- lapply(1:n.perm, function(p) one.perm(p))	
					
				}
			#also make variables to hold the permutation results
			good.results.perm <- which(unlist(lapply(perm.results, function(x) length(x$model.effects))) > 0)
			perm.results <- perm.results[good.results.perm]			
			all.model.effects.perm <- matrix(unlist(lapply(perm.results, function(x) x$model.effects)), nrow = length(good.results.perm), byrow = TRUE)
			all.model.se.perm <- matrix(unlist(lapply(perm.results, function(x) x$model.se)), nrow = length(good.results.perm), byrow = TRUE)
			all.model.cov.perm <- matrix(unlist(lapply(perm.results, function(x) x$model.cov)), nrow = length(good.results.perm), byrow = TRUE)
			marker.pairs.used.perm <- matrix(unlist(lapply(perm.results, function(x) x$pair.used)), nrow = length(good.results.perm), byrow = TRUE)
			colnames(all.model.effects.perm) <- colnames(all.model.se.perm) <- column.names			
			rm(perm.results)
		      
      
			perm.paired.marker.num <- matrix(apply(marker.pairs.used.perm, 2, function(x) data.obj$marker.num[match(x, data.obj$marker.names)]), ncol = 2, byrow = FALSE)
	      	colnames(perm.paired.marker.num) <- c("marker1", "marker2")
			final.effects.table.perm <- cbind(perm.paired.marker.num, all.model.effects.perm)
			final.se.table.perm <- cbind(perm.paired.marker.num, all.model.se.perm)	
			
			final.cov.table.perm <- all.model.cov.perm
			phenotype.perm.results <- list(final.effects.table.perm, final.se.table.perm, final.cov.table.perm)
			names(phenotype.perm.results) <- c("pairscan.effects.perm", "pairscan.se.perm", "model.covariance.perm")
      		}else{
      			phenotype.perm.results <- NULL
      			}
	      	
	      	
	      	final.results <- list(phenotype.results, phenotype.perm.results)
	      	names(final.results) <- c("pairscan.results", "pairscan.perm")

			if(verbose){
				cat("\n") #make sure the prompt is on the next line at the end of everything
				}


	      	return(final.results)


	}
