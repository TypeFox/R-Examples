get.pairs.for.pairscan <-
function(geno, min.per.genotype = NULL, max.pair.cor = NULL, verbose = FALSE, n.cores = 2){

	p = NULL #for appeasing R CMD check

	if(is.null(min.per.genotype) && is.null(max.pair.cor)){
		stop("Either min.per.genotype or max.pair.cor should be set. Type ?get.pairs.for.pairscan for help.")
		}
	
	if(!is.null(min.per.genotype) && !is.null(max.pair.cor)){
		stop("Only one of min.per.genotype or max.pair.cor can be specified. Type ?get.pairs.for.pairscan for help.")
		}
	
	if(verbose){
		cat("\nChecking marker pairs for genotype representation...\n")
		}
	
	
	if(!is.null(min.per.genotype)){
		thresh.param <- min.per.genotype
		check.linkage <- function(m1,m2,thresh.param){
			t <- cbind.data.frame(as.factor(m1),as.factor(m2))
			colnames(t) <- c("m1","m2")
			reps <- table(t$m1,t$m2)
			too.few <- which(reps < min.per.genotype)
			if(length(too.few) >= 1) {
				return(FALSE) #pair failed check
				}else{
				return(TRUE) #pair passed check
				}
			}
		}else{
		thresh.param <- max.pair.cor
		check.linkage <- function(m1,m2,thresh.param){
			pair.cor <- try(cor(m1, m2, use = "complete"), silent = TRUE)
			if(class(pair.cor) == "try-error" || pair.cor > max.pair.cor) {
				return(FALSE) #pair failed check
				}else{
				return(TRUE) #pair passed check
				}
			}
		}


	
	all.pairs <- pair.matrix(1:dim(geno)[2])
	all.pair.names <- pair.matrix(colnames(geno))
	
	check.one.pair <- function(pair){
		pass.checks <- check.linkage(m1 = geno[,pair[1]], m2 = geno[,pair[2]], thresh.param = thresh.param)
		return(pass.checks)
		}
	
	cl <- makeCluster(n.cores)
	registerDoParallel(cl)
	good.pairs <- foreach(p = t(all.pairs), .combine = "c") %dopar% {
		check.one.pair(p)
		}
	stopCluster(cl)

		
	pairs.mat <- all.pair.names[which(good.pairs),,drop = FALSE]
	colnames(pairs.mat) <- c("marker1", "marker2")
	rownames(pairs.mat) <- NULL
	
	if(verbose){
		cat(dim(pairs.mat)[1], "pairs will be tested.\n")
		}
	
	return(pairs.mat)
	
}
