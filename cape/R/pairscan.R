pairscan <-
function(data.obj, geno.obj = NULL, covar = NULL, scan.what = c("eigentraits", "raw.traits"), total.perm = NULL, min.per.genotype = NULL, max.pair.cor = NULL, n.top.markers = NULL, use.kinship = FALSE, kin.full.geno = TRUE, sample.kinship = TRUE, num.kin.samples = 100, n.per.sample = 100, verbose = FALSE, num.pairs.limit = 1e6, num.perm.limit = 1e7, overwrite.alert = TRUE, n.cores = NULL) {

	covar.type = NULL #for appeasing R CMD check


	if(overwrite.alert){
		choice <- readline(prompt = "Please make sure you assign the output of this function to a pairscan.obj, and NOT the data.obj. It will overwrite the data.obj.\nDo you want to continue (y/n) ")
		if(choice == "n"){stop()}
		}


	pairscan.obj <- list()

	if(is.null(total.perm)){
		stop("The number of permutations must be specified.")
		}
	
	pheno <- get.pheno(data.obj, scan.what)	
	num.pheno <- dim(pheno)[2]
	pheno.names <- colnames(pheno)

	covar.info <- get.covar(data.obj, covar)
	for(i in 1:length(covar.info)){
		assign(names(covar.info)[i], covar.info[[i]])
		}
	#find the phenotypic covariates. These will
	#be tested separately, and not as part of a
	#chromosome
	num.covar <- length(covar.type)
	p.covar.locale <- which(covar.type == "p")
	num.p.covar <- length(p.covar.locale)

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
	add.full.kin <- as.logical(num.p.covar > 0)
	

	if(is.null(data.obj$geno.for.pairscan)){
		stop("select.markers.for.pairscan() must be run before pairscan()")
		}
		
	
	#add the covariates (geno and pheno) 
	#to the genotype matrix so that we 
	#test all pairs
	geno <- get.geno.with.covar(data.obj, geno.obj, g.covar = TRUE, p.covar = TRUE, for.pairscan = TRUE)	
	num.markers <- dim(geno)[2]
	
	#fill in a matrix to index the marker pairs
	pared.marker.mat <- get.pairs.for.pairscan(geno, min.per.genotype = min.per.genotype, max.pair.cor = max.pair.cor, verbose = verbose, n.cores = n.cores)

	num.pairs <- dim(pared.marker.mat)[1]
	
	if(num.pairs == 0){
		stop("There are no pairs to test. Try lowering min.per.genotype or raising max.pair.cor.")
		}

	if(!is.null(num.pairs.limit) && num.pairs > num.pairs.limit){
		cat("\nThe number of pairs (",num.pairs,") exceeds ", num.pairs.limit, ".\n", sep = "")
		go.on <- readline(prompt = "Do you want to continue (y/n)?\n")
		if(length(grep("n", go.on))){
			message("Stopping pairwise scan...\n")
			return(pairscan.obj)
		}else{
			cat("Continuing pairwise scan...\n")
		}
	}


	#run one.pairscan for each phenotype with results in scanone.result
	if(!use.kinship){
		pairscan.results <- pairscan.noKin(data.obj = data.obj, pheno.mat = pheno, geno.mat = geno, covar.table = covar.table, paired.markers = pared.marker.mat, n.perm = total.perm, verbose = verbose, n.cores = n.cores)
		}else{
		pairscan.results <- pairscan.kin(data.obj = data.obj, geno.obj = geno.obj, covar = covar, scan.what = scan.what, marker.pairs = pared.marker.mat, kin.full.geno = kin.full.geno, sample.kinship = sample.kinship, num.kin.samples = num.kin.samples, n.per.sample = n.per.sample, verbose = verbose, n.cores = n.cores)
		}	
		
	pairscan.obj$pairscan.results <- pairscan.results	
	
	#calculate the number of top markers to use based
	#on the thresholding of the singlescan
	if(is.null(n.top.markers)){
		n.top.markers <- dim(geno)[2]
		}

	pairscan.perm <- pairscan.null(data.obj, geno.mat = geno, covar = covar, scan.what = scan.what, total.perm = total.perm, n.top.markers = n.top.markers, verbose = verbose, min.per.genotype = min.per.genotype, max.pair.cor = max.pair.cor, n.cores = n.cores)


	pairscan.obj$pairscan.perm <- pairscan.perm$pairscan.perm #add the results to the data object
	pairscan.obj$pairs.tested.perm <- pairscan.perm$pairs.tested.perm #add the results to the data object


	return(pairscan.obj) #and return it
	
}
