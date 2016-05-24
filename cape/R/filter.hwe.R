filter.hwe <-
function(data.obj, geno.obj, p.thresh = 1e-6, run.parallel = TRUE, n.cores = 2){
	# library(HardyWeinberg)
	
	m = NULL #for appeasing R CMD check

	geno.mat <- get.geno(data.obj, geno.obj)
	
	removed.markers <- NULL
	
	get.hwe <- function(x){
		genotypes <- rep(NA, 3)
		names(genotypes) <- c("AA", "AB", "BB")
		genotypes[1] <- length(which(x == 0))
		genotypes[2] <- length(which(x == 0.5))
		genotypes[3] <- length(which(x == 1))
		hwe.p <- HWExact(genotypes)$pval
		return(hwe.p)
		}
	
	
	if(run.parallel){
		# par.start <- proc.time()
		cl <- makeCluster(n.cores)
		registerDoParallel(cl)
		all.hwe <- foreach(m = geno.mat, .combine = "c") %dopar% {
			get.hwe(m)
			}
		stopCluster(cl)
		
		# par.end <- proc.time()	
		}else{
		# ser.start <- proc.time()
		all.hwe <- apply(geno.mat, 2, get.hwe)
		# ser.end <- proc.time()
		}
			
	# par.time <- par.end - par.start
	# ser.time <- ser.end - ser.start
			
	too.low <- which(all.hwe < p.thresh)
	if(length(too.low) > 0){
		geno.mat <- geno.mat[,-too.low]
		removed.markers <- c(removed.markers, names(too.low))
		}
	
	if(length(removed.markers) > 0){
		cat("Removed", length(removed.markers), "markers due to deviation from HWE.\n")
		write.table(cbind(sort(removed.markers), data.obj$marker.names[sort(removed.markers)]), "markers.removed.for.HWE.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
		data.obj <- remove.markers(data.obj, removed.markers)
		}
	
	return(data.obj)
	}
