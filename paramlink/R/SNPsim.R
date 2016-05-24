.SNPsim <- function(x, N=1, partialmarker=NULL, available=x$available, afreq=c(0.5, 0.5), loop_breakers=NULL, unique=FALSE, seed=NULL, verbose=TRUE) {
	if (is.null(x$model)) {
		if (all(x$pedigree[,'AFF']==1)) {x=setModel(x,1); if(verbose) cat("Unaffected pedigree. Simulating autosomal markers.\n")}
		else stop("No model set. Use setModel().")
	}
	
	if (any(!is.numeric(N), length(N)>1, N%%1 != 0)) stop("N must be a positive integer.")
	if (length(afreq)>2) stop("Sorry - only diallelic markers allowed in SNP simulation.")
	nInd = x$nInd; chrom = x$model$chrom; 
	
	.TRzero <- .TRmatr(0, chrom); .TRhalf <- .TRmatr(0.5, chrom)
	
	if (is.null(partialmarker)) partialmarker = matrix(0, nrow=nInd, ncol=2) else if(ncol(partialmarker)!=2 || nrow(partialmarker)!=nInd) stop("Wrong dimensions of marker matrix.")
	x = setMarkers(x, partialmarker, missing=0)
	if(length(available) == 0) available = x$orig.ids
	x$available = available  #this must be reset because of setMarkers()
	
	m = x$markerdata[[1]]
	if(length(attr(m, "alleles")[[1]])>2) stop("Partial marker is not diallelic. Only diallelic markers allowed.")
	if(verbose && any(m!=0)) { 
		cat("Simulating markers conditional on existing genotypes:\n")
		print(data.frame(ID=x$orig.ids, GENO=.prettyMarkers(m, missing="-", singleCol=TRUE, sex=x$pedigree[, 'SEX'])), row.names=FALSE)
	}

	if (loops <- x$hasLoops)	{		
		if(is.null(loop_breakers)) 	 stop("The pedigree has loops. Please indicate loop breakers.")
		x = breakLoops(x, loop_breakers)
		m = x$markerdata[[1]]	
		nInd = x$nInd
	}

	zgeno = .diallel2geno(m); preexisting = zgeno!=0
	# create initial marker matrix: one column per marker (single-numerically coded)
	markers = matrix(rep(zgeno, N), ncol=N)
	
	avail = x$orig.ids %in% x$available
	if (!is.null(seed)) set.seed(seed)
	
	# simulate only indivs i) with no given genotype, ii) that are available or have available descendants. Duplicated individuals are not simulated.
	sim_indivs = seq_len(nInd)[sapply(seq_len(nInd), function(i) zgeno[i]==0 && (avail[i] || any(avail[descendants(x, x$orig.ids[i])])))]
	if(length(sim_indivs)==0) stop("Something is wrong: No individuals available for simulation.")
	if(verbose) cat("Simulating genotypes for the following individuals:", .prettycat(x$orig.ids[avail & !preexisting], "and"), "\n")
	
	genoprobs <- function(x, partialgeno, id, values) {  
			#values are 1:3 for autosomal models, and either 1:2 (males) or 1:3 (females) for X-linked models
			#outputs vector of length |values| with genotype probs for indiv id given pedigree og partial genotype information
			probs <- sapply(values, function(g) { partialgeno[id]<-g;   .likelihoodSNP(x, afreq=afreq, singleNum.geno=partialgeno, TR.MATR=.TRzero) } )
			if (sum(probs)==0) { print(cbind(x$pedigree, partialgeno)); stop("\nIndividual ",id,": All genotype probabilities zero. Mendelian error?") }
			probs
	}
	
	switch(chrom,	
	AUTOSOMAL = {
		# simulation order: founders first (speeds up the likelihoods). 
		sim_indivs <- sim_indivs[order(sim_indivs %in% x$nonfounders)]
		
		# pre-calculate probabilities for the first 'init' individuals (big time saver!)
		init <- which.min(3^seq_along(sim_indivs) + 3*N*(length(sim_indivs)-seq_along(sim_indivs)))  #optimal 'init' minimizes the number of likelihood() calls
		initg <- t.default(fast.grid( rep(list(1:3), init ) ))
		initp <- apply(initg, 2, function(g) { 
			zgeno[ sim_indivs[1:init] ] <- g;   
			.likelihoodSNP(x, afreq=afreq, singleNum.geno=zgeno, TR.MATR=.TRzero) } )
		if (identical(sum(initp), 0)) stop("All genotype probabilities zero. Wrong model?")
		
		# pre-fill the rows of the 'init' individuals (i.e. sim_indivs[1:init]) 
		markers[sim_indivs[1:init], ] <- initg[ , suppressWarnings(sample.int( ncol(initg), size=N, replace=TRUE, prob=initp )) ]

		# do the rest of the individuals (only those present in sim_indivs)
		for (i in sim_indivs[-(1:init)])
			markers[i, ] <- apply(markers, 2, function(partgeno) sample.int(3, size=1, prob=genoprobs(x, partgeno, id=i, values=1:3)) )
	},
	X = {
		ped = x$pedigree
		males = sim_indivs[ped[,'SEX'][sim_indivs]==1]; females = sim_indivs[ped[,'SEX'][sim_indivs]==2]
		males = males[order(males %in% x$nonfounders)]; females = females[order(females %in% x$nonfounders)]  #quicker with this?
		n_males=length(males); n_females=length(females)
		
		# find optimal 'init' values for males/females
		calls = outer(0:n_males, 0:n_females, function(ma, fe) 2^ma * 3^fe + 2*(n_males-ma)*N + 3*(n_females-fe)*N) # = number of times likelihood() is called.
		calls.min = arrayInd(which.min(calls), dim(calls)) 
		init_m = calls.min[1]-1; init_f = calls.min[2]-1
		init_indivs = c(males[seq_len(init_m)], females[seq_len(init_f)])
		# pre-calculate probabilities for the first 'init' individuals
		initg <- t.default(fast.grid( list(1:2, 1:3)[rep(1:2, c(init_m, init_f))] ))
		initp <- apply(initg, 2, function(g) { zgeno[ init_indivs ] <- g;   .likelihoodSNP(x, afreq=afreq, singleNum.geno=zgeno, TR.MATR=.TRzero) } )

		# pre-fill the rows of the 'init_indivs' 
		markers[init_indivs, ] <- initg[ , suppressWarnings(sample.int( ncol(initg), size=N, replace=TRUE, prob=initp)) ]

		# do the rest of the males (only those present in sim_indivs)
		for (i in setdiff(males, males[seq_len(init_m)])) 
			markers[i, ] <- apply(markers, 2, function(partgeno) sample.int(2, size=1, prob=genoprobs(x, partgeno, id=i, values=1:2)) )
		
		for (i in setdiff(females, females[seq_len(init_f)])) 
			markers[i, ] <- apply(markers, 2, function(partgeno) sample.int(3, size=1, prob=genoprobs(x, partgeno, id=i, values=1:3)) ) }
	)
	
	markers[!avail & !preexisting, ] <- 0
	
	if (unique) markers = unique(markers, MARGIN=2)
	markers2 = .geno2diallel(markers)
	mlist = lapply(2*seq_len(ncol(markers)), function(i) .createMarkerObject(markers2[, (i-1):i], alleles=1:2, afreq = afreq, missing = 0))
	x = setMarkers(x, structure(mlist, class = "markerdata"))	
	if(loops) x = tieLoops(x)
	x
}
