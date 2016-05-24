impute.missing.geno <-
function(data.obj, geno.obj = NULL, impute.full.genome = FALSE, k = 10, ind.missing.thresh = 0, marker.missing.thresh = 0, prioritize = c("ind", "marker", "fewer"), max.region.size = NULL, min.region.size = NULL, run.parallel = TRUE, verbose = FALSE, n.cores = 2){
	
		
	p.choice <- grep("f", prioritize)

	if(length(p.choice) > 0){prioritize <- "fewer"}
	
	if(!impute.full.genome){
		new.geno <- get.geno(data.obj, geno.obj)
		chr <- data.obj$chromosome
		u_chr <- unique(chr)
		}else{
		new.geno <- geno.obj$geno	
		if(is.null(new.geno)){
			stop("I could not find a genotype matrix in the genotype object. Please check or set impute.full.genome to FALSE")
			}
		chr <- geno.obj$chromosome
		u_chr <- unique(chr)
		ind.used <- rownames(data.obj$pheno)
		ind.used.locale <- match(ind.used, rownames(new.geno))
		new.geno <- new.geno[ind.used.locale,]
		}

	

		chr.lengths <- apply(matrix(u_chr, nrow = 1), 2, function(x) length(which(chr == x)))
		
		#go through the chromosome at the window size specified and impute genotypes
		if(is.null(max.region.size) || max.region.size > max(chr.lengths)){
			max.region <- max(chr.lengths)
			min.region <- min(chr.lengths)
			}else{
			max.region <- max.region.size
			min.region <- min.region.size
			if(is.null(min.region)){
				min.region <- 10			
				}
			}
		if(max.region > 1000){
			choice <- 0
			while(choice != "y" && choice != "n"){
			choice <- readline(prompt = "At least one chromosome has more than 1000 markers.\n\tCalculation of large distance matrices may crash R.\n\tSet max.region.size to specify maximum number of markers in imputed region.\n\tWould you like to continue with the current settings? (y/n)")
			if(choice == "n"){stop()}
			}	
		}


	marker.names <- colnames(new.geno)
	ind.names <- rownames(new.geno)
	
	num.sections <- length(seq(1, dim(new.geno)[2], max.region))
	geno.chunks <- vector(mode = "list", length = num.sections)

	if(verbose){cat("Dividing genome into chunks.\n")}
	chunk.pos <- 1
	for(ch in u_chr){
		if(verbose){report.progress(which(u_chr == ch), length(u_chr))}
		#pull out all the genotypes for the chromosome
		chr.locale <- which(chr == ch)
		if(length(chr.locale) < max.region){
			geno.chunks[[chunk.pos]] <- new.geno[,chr.locale,drop=FALSE]
			chunk.pos <- chunk.pos + 1
			}else{
			marker.seq <- seq(1, length(chr.locale), max.region)
			if(tail(marker.seq, 1) < length(chr.locale)){
				if(length(chr.locale)- tail(marker.seq,1) >= min.region){
					marker.seq <- c(marker.seq, length(chr.locale))
					}else{
					marker.seq[length(marker.seq)] <- length(chr.locale)
					}
				}
			for(m in 1:(length(marker.seq)-1)){
				if(m == length(marker.seq)-1){
					last.pos <- chr.locale[(marker.seq[m+1])]
					}else{
					last.pos <- chr.locale[(marker.seq[m+1]-1)]	
					}
				geno.chunks[[chunk.pos]] <- new.geno[,chr.locale[marker.seq[m]]:last.pos]
				chunk.pos <- chunk.pos + 1
				}
			}
		}
	if(verbose){cat("\n")}
		
	impute.section <- function(sec.geno){
		if(dim(sec.geno)[2] == 1){return(sec.geno)}
		ind.missing.geno <- which(apply(sec.geno, 1, function(x) length(which(is.na(x)))) > 0)
			if(length(ind.missing.geno) > 0){
				ind.dist <- as.matrix(dist(sec.geno))
				for(id in 1:length(ind.missing.geno)){
					ind <- sec.geno[ind.missing.geno[id],,drop=FALSE]
					neighbors <- sort(ind.dist[ind.missing.geno[id],])
					if(length(neighbors) < 2){
						next()
						}
					nearest.neighbors <- as.numeric(names(neighbors[2:min(c((k+1), length(neighbors)))]))
					neighbor.weights <- neighbors[nearest.neighbors]/sum(neighbors[nearest.neighbors], na.rm = TRUE)
					missing.geno <- which(is.na(ind))
					missing.geno.ind <- sec.geno[as.character(nearest.neighbors),missing.geno,drop=FALSE]
					weight.mat <- matrix(neighbor.weights, ncol = ncol(missing.geno.ind), nrow = nrow(missing.geno.ind), byrow = FALSE)
					missing.geno.weighted <- missing.geno.ind*weight.mat
					ind[missing.geno] <- apply(missing.geno.weighted, 2, function(x) sum(x, na.rm = TRUE))
					sec.geno[ind.missing.geno[id],] <- ind
					} #end looping through individuals with missing genotypes in this section
				} #end case for if there are individuals with missing genotypes
			return(sec.geno)
		} #end impute.section
	
				
		if(run.parallel){		
			if(verbose){cat("Imputing missing genotypes...\n")}
			cl <- makeCluster(n.cores)
			registerDoParallel(cl)
			imputed.geno <- foreach(m = geno.chunks) %dopar% {
				impute.section(m)
				}
				stopCluster(cl)
			}else{
			# if(verbose){cat("Imputing missing genotypes...\n")
				# imputed.geno <- lapply_pb(geno.chunks, impute.section)
				# }else{
				imputed.geno <- lapply(geno.chunks, impute.section)
				# }
			}
	
		if(verbose){cat("Rebuilding genotype matrix...\n")}
		#now unlist the result into the new genotype matrix
		new.col <- sum(unlist(lapply(imputed.geno, function(x) dim(x)[2])))
		imp.geno <- matrix(unlist(imputed.geno), ncol = new.col, byrow = FALSE)
		colnames(imp.geno) <- marker.names; rownames(imp.geno) <- ind.names
		
	#========================================================================
	# internal functions
	#========================================================================
	assess.missing <- function(geno.check){
		ind.missing.percent <- as.vector(apply(geno.check, 1, function(x) length(which(is.na(x)))))/dim(geno.check)[2]*100
		marker.missing.percent <- as.vector(apply(geno.check, 2, function(x) length(which(is.na(x)))))/dim(geno.check)[1]*100
		
		ind.missing.lots <- which(ind.missing.percent > ind.missing.thresh)
		marker.missing.lots <- which(marker.missing.percent > marker.missing.thresh)
		
		results <- list(ind.missing.lots, marker.missing.lots)
		names(results) <- c("ind.missing.lots", "marker.missing.lots")
		return(results)
		}	
	
	remove.ind.int <- function(data.obj, geno, missing.ind){
		ind <- missing.ind$ind.missing.lots
		ind.id <- rownames(geno)[ind]
		ind.locale <- match(ind.id, rownames(data.obj$pheno))
		if(length(ind) > 0){
			cat(paste("Removing ", length(ind), " individual(s) with more than ", ind.missing.thresh, "% missing data.\n", sep = ""))
			#remove individuals from the imputed genotype matrix
			#and from the phenotype matrix in data.obj
			geno <- geno[-ind,,drop=FALSE]
			data.obj$pheno <- data.obj$pheno[-ind.locale,,drop=FALSE]
			if(!is.null(data.obj$p.covar.table)){
				data.obj$p.covar.table <- data.obj$p.covar.table[-ind.locale,,drop=FALSE]
				}
			if(!is.null(data.obj$g.covar.table)){
				data.obj$g.covar.table <- data.obj$g.covar.table[-ind.locale,,drop=FALSE]
				}
			if(!is.null(data.obj$raw.pheno)){
				data.obj$raw.pheno <- data.obj$raw.pheno[-ind.locale,,drop=FALSE]
				}
			if(!is.null(data.obj$ET)){
				warning("get.eigentraits needs to be re-run because individuals were removed.\n")	
				}
			}
		results.list <- vector(mode = "list", length = 2)
		names(results.list) <- c("data.obj", "geno.obj")
		results.list[[1]] <- data.obj
		results.list[[2]] <- geno
		return(results.list)
		}


	#This function removes markers from the imputed genotype
	#matrix, updates the meta information in the data.obj
	#if there is a genotype matrix in the data.obj, it is
	#replaced with an imputed version
	remove.marker <- function(data.obj, geno, missing.ind){
		marker.missing.lots <- missing.ind[[2]]
		marker.id <- colnames(geno)[marker.missing.lots]
		
		if(length(marker.missing.lots) > 0){
			cat(paste("Removing ", length(marker.missing.lots), " markers with more than ", marker.missing.thresh, "% missing data.\n", sep = ""))
			geno <- geno[,-marker.missing.lots,drop=FALSE]
			}

		#if any of these markers are in the data.obj, take them out
		marker.d.ind <- which(data.obj$marker.num %in% marker.id)
		
		#take out the markers from the data.obj meta data
		if(length(marker.d.ind) > 0){
			data.obj$chromosome <- data.obj$chromosome[-marker.d.ind]
			data.obj$marker.names <- data.obj$marker.names[-marker.d.ind]
			data.obj$marker.num <- data.obj$marker.num[-marker.d.ind]
			data.obj$marker.location <- data.obj$marker.location[-marker.d.ind]
			}
		
		#if there is a genotype matrix in the data.obj, update it with
		#the imputed genotype matrix
		if(!is.null(data.obj$geno)){
			data.obj <- update.data.obj(data.obj, geno)
			}
			
		results.list <- vector(mode = "list", length = 2)
		names(results.list) <- c("data.obj", "geno.obj")
		results.list[[1]] <- data.obj
		results.list[[2]] <- geno
		return(results.list)
		}
		
		#This function is used to replace the genotype matrix
		#in the geno.obj with the imputed genotypes
		update.geno.obj <- function(geno.obj, geno){			
			#locate the imputed marker names in the geno.obj
			marker.names <- 
			imputed.marker.locale <- match(colnames(geno), geno.obj$marker.names)
			imputed.ind.locale <- match(rownames(geno), rownames(geno.obj$geno))
			geno.obj$geno[imputed.ind.locale,imputed.marker.locale] <- imp.geno
			geno.obj$marker.names <- geno.obj$marker.names[imputed.marker.locale]
			geno.obj$chromosome <- geno.obj$chromosome[imputed.marker.locale]
			geno.obj$marker.location <- geno.obj$marker.location[imputed.marker.locale]
			geno.obj$marker.num <- geno.obj$marker.num[imputed.marker.locale]
			return(geno.obj)
			}
		
		#This function is used to put the imputed genotypes into
		#the data.obj if that's where the genotype data are stored
		update.data.obj <- function(data.obj, geno){
			marker.locale <- get.marker.idx(data.obj, colnames(geno))
			ind.locale <- match(rownames(data.obj$pheno), rownames(geno))
			data.obj$geno <- geno[ind.locale,marker.locale]
			colnames(data.obj$geno) <- data.obj$marker.num
			rownames(data.obj$geno) <- rownames(data.obj$pheno)
			return(data.obj)
			}
		#========================================================================
		# end internal functions
		#========================================================================
		
	if(verbose){cat("Assessing missing values...\n")}
	num.missing.geno <- length(which(is.na(imp.geno)))
	
	if(num.missing.geno > 0){
		if(prioritize == "ind"){
			test <- assess.missing(imp.geno)
			if(length(test$ind.missing.lots) > 0){
				results <- remove.ind.int(data.obj, imp.geno, test)
				}
			if(length(test$marker.missing.lots) > 0){
				results2 <- remove.marker(results[[1]], results[[2]], test)
				}
			data.obj <- results2[[1]]
			imp.geno <- results2[[2]]
			}
	
		if(prioritize == "marker"){
			test <- assess.missing(imp.geno)
			results <- remove.marker(data.obj, geno = imp.geno, missing.ind = test)
			results2 <- remove.ind.int(results[[1]], results[[2]], test)
			data.obj <- results2[[1]]
			imp.geno <- results2[[2]]
			}
			
		if(prioritize == "fewer"){
			test <- assess.missing(imp.geno)
			if(length(test$ind.missing.lots) < length(test$marker.missing.lots)){
				results <- remove.ind.int(data.obj, geno = imp.geno, missing.ind = test)
				results2 <- remove.marker(data.obj = results[[1]], geno = results[[2]], missing.ind = test)
				data.obj <- results2[[1]]
				imp.geno <- results2[[2]]
				}else{
				results <- remove.marker(data.obj, imp.geno, test)
				results2 <- remove.ind.int(results[[1]], results[[2]], test)
				data.obj <- results2[[1]]
				imp.geno <- results2[[2]]
				}
			}
		
		
		#the data.obj is now finalized. If a geno.obj is provided
		#it needs to be updated with the imputed genotypes
		if(!is.null(data.obj$geno)){
			#if we did not impute the full genome, 
			#oly return the data object
			if(impute.full.genome == FALSE){
				geno.obj <- NULL
				results.list <- vector(mode = "list", length = 2)
				names(results.list) <- c("data.obj", "geno.obj")
				results.list[[1]] <- data.obj
				results.list[[2]] <- geno.obj
				return(results.list)
				# return(data.obj)
				}else{
				#if we did, update the geno.obj and return
				#both the data.obj and the geno.obj
				geno.obj <- update.geno.obj(geno.obj, imp.geno)
				results.list <- vector(mode = "list", length = 2)
				names(results.list) <- c("data.obj", "geno.obj")
				results.list[[1]] <- data.obj
				results.list[[2]] <- geno.obj
				return(results.list)
				}
			}else{
			geno.obj <- update.geno.obj(geno.obj, imp.geno)
			results.list <- vector(mode = "list", length = 2)
			names(results.list) <- c("data.obj", "geno.obj")
			results.list[[1]] <- data.obj
			results.list[[2]] <- geno.obj
			return(results.list)
			}
		}else{ 
			#if there are no missing values after imputation, 
			#return the data.obj, and the geno.obj with the
			#imputed genotypes
			if(!is.null(data.obj$geno)){
				if(impute.full.genome == FALSE){
					#if we are only imputing the genotypes
					#in the data.obj, just put them back
					data.obj$geno <- imp.geno
					geno.obj <- NULL
					results.list <- vector(mode = "list", length = 2)
					names(results.list) <- c("data.obj", "geno.obj")
					results.list[[1]] <- data.obj
					results.list[[2]] <- geno.obj
					return(results.list)
					}else{
					#if we are imputing the full genome
					#put only the correct markers back in the data.obj
					data.obj <- update.data.obj(data.obj, imp.geno)
					geno.obj <- update.geno.obj(geno.obj, imp.geno)
					results.list <- vector(mode = "list", length = 2)
					names(results.list) <- c("data.obj", "geno.obj")
					results.list[[1]] <- data.obj
					results.list[[2]] <- geno.obj
					return(results.list)
					}
				} #end case for if there is a genotype matrix in the data.obj
				#if there is no genotype matrix in the data.obj, we update
				#the geno.obj and return both the data.obj and the geno.obj
				geno.obj <- update.geno.obj(geno.obj, geno = imp.geno)
				results.list <- vector(mode = "list", length = 2)
				names(results.list) <- c("data.obj", "geno.obj")
				results.list[[1]] <- data.obj
				results.list[[2]] <- geno.obj
				return(results.list)
				return(list(data.obj, geno.obj))
				}

}
