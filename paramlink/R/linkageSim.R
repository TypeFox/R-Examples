linkageSim <- function(x, N=1, available=x$available, afreq=NULL, partialmarker=NULL, 
                   loop_breakers=NULL, unique=FALSE, seed=NULL, verbose=TRUE) {
    if (any(!is.numeric(N), length(N)>1, N%%1 != 0)) stop("N must be a positive integer.")
    if (is.null(x$model)) stop("No model set. Use setModel().")
    if ((loops <- x$hasLoops) && is.null(loop_breakers))  stop("The pedigree has loops. Please indicate loop breakers.")
    if (!is.null(partialmarker)) {
      stop("Not implemented yet")
        if (!inherits(partialmarker, "marker")) stop("Argument 'partialmarker' must be a 'marker' object (or NULL).")
        else if (nrow(partialmarker)!=x$nInd) stop("Partial marker does not fit the pedigree.")
        if(!is.null(afreq)) stop("When 'partialmarker' is non-NULL, 'afreq' must be NULL.")
        if(length(mendelianCheck(setMarkers(x, partialmarker), verbose=F)) > 0) stop("Mendelian error in the given partial marker.")
    }
    else {
        if(is.null(afreq)) stop("Please specify marker allele frequencies.")
        partialmarker = marker(x, alleles=seq_along(afreq), afreq=afreq)
    }
    if (!is.null(seed)) 
        set.seed(seed)
    starttime = proc.time()
    
    likel_counter = 0
    m = partialmarker
    afreq = attr(m, 'afreq')
    chrom = x$model$chrom; nall = length(afreq);
    if(length(available) == 0) 
        available = x$orig.ids
    sim_origs = available
   
    if (nall > 4)
        stop("Sorry - simulation of markers with more than 4 alleles is not implemented yet.")
    
    if(verbose && any(m!=0)) { #not implemented
        cat("Simulating markers conditional on existing genotypes:\n")
        print(data.frame(ID=x$orig.ids, GENO=.prettyMarkers(m, missing="-", singleCol=TRUE, sex=x$pedigree[, 'SEX'])), row.names=FALSE)
    }

    if(loops) {        
        lb = loop_breakers
        x = breakLoops(setMarkers(x, m), lb)
        #if(verbose) cat(ifelse(length(lb)==1, "Breaking loop at individual ", "\nBreaking loops at individuals "), 
        #                .prettycat(lb, "and"), "\n", sep="")
        m = x$markerdata[[1]]    
        sim_origs = unique.default(c(sim_origs, lb))
    }
   
    initialCalc = .initialCalc(x, afreq, chrom)

    # simulation order: founders first, and those with many possible haplotypes (suggested by quick tests...)
    sim_indivs = .internalID(x, sim_origs)
    loop_int = .internalID(x, x$loop_breakers[, 1])
    
    sim_indivs = sim_indivs[order(!sim_indivs %in% loop_int, sim_indivs %in% x$nonfounders)]#, initz)]

    if(length(sim_indivs)==0) stop("Something is wrong: No individuals available for simulation.")
    if(verbose) {
      cat("Target individuals:", .prettycat(available, "and"), "\n")
      if(length(rest <- setdiff(sim_origs, available)) > 0) 
         cat("Additional individuals simulated internally (increasing speed):", .prettycat(rest, "and"), "\n")
      
    }
    .TRzero = .TRmatrNEW(0, nall, chrom)
    
    gt_L = nall*(nall+1)/2 #number of genotype codes (= # unordered genotypes)
    gt_values = seq_len(gt_L)
   
    # create initial marker matrix: one column per marker (single-numerically coded)
    zgeno = .diallel2genoNEW(m, nall)
    markers = matrix(rep(zgeno, N), ncol=N)
    
    # pre-calculate probabilities for the first 'init' individuals (big time saver!)
    switch(chrom,    
    AUTOSOMAL = {
        calls = gt_L^seq_along(sim_indivs) + gt_L*N*(length(sim_indivs)-seq_along(sim_indivs))
        init = which.min(calls)  #optimal 'init' minimizes the number of likelihood() calls
        if(init==0) stop("init = 0")
        init_int = sim_indivs[1:init]
        cost_int = sim_indivs[-(1:init)]
        initgrid = fast.grid( rep(list(gt_values), init ) )
    },
    X = {
        SEX = x$pedigree[,'SEX']
        males_int = sim_indivs[SEX[sim_indivs]==1]; females_int = sim_indivs[SEX[sim_indivs]==2]
        n_males=length(males_int); n_females=length(females_int)
        loop_int = .internalID(x, loop_breakers); loop_m = sum(SEX[loop_int]==1); loop_f = sum(SEX[loop_int]==2)
      
        # find optimal 'init' values for males/females
        gtM = nall;  gt_values_m = 1:nall; gtF = gt_L
        calls = outer(loop_m:n_males, loop_f:n_females, function(ma, fe)      # = number of times likelihood() is called.
               gtM^ma * gtF^fe + gtM*(n_males-ma)*N + gtF*(n_females-fe)*N) 
        
        calls.min = arrayInd(which.min(calls), dim(calls)) 
        init_m = calls.min[1] - 1 + loop_m; init_f = calls.min[2] - 1 + loop_f
        init_int = c(males_int[seq_len(init_m)], females_int[seq_len(init_f)])
        cost_int = setdiff(sim_indivs, init_int)
        initgrid = fast.grid( list(gt_values_m, gt_values)[rep(1:2, c(init_m, init_f))] )
    })
   
    if(verbose) cat("\nSimulation order:\n   Precomputing joint probabilities:", .prettycat(x$orig.ids[init_int], 'and'),
                        "\n   Brute force:", .prettycat(x$orig.ids[cost_int], 'and'), "\n")
     
    # pre-fill the rows of the 'init' individuals (i.e. sim_indivs[1:init]) 
    initp = apply(initgrid, 1, function(g) { 
        zgeno[init_int] = g
        likelihood_LINKAGE(x, afreq=afreq, singleNum.geno=zgeno, initialCalc=initialCalc, TR.MATR=.TRzero) 
    })
    likel_counter = likel_counter + length(initp)
    if (identical(sum(initp), 0)) stop("All genotype probabilities zero. Wrong model?")
   
    markers[init_int, ] = t.default(initgrid[ suppressWarnings(sample.int( nrow(initgrid), size=N, replace=TRUE, prob=initp)), ])
    if(verbose) cat("\nPrecomputing finished. Time used:", (proc.time()-starttime)[3], "seconds.\n")
   

    genoprobs <- function(x, partial, id, values) {  
        #values are 1:gt_L for autosomal models, and either 1:nall (males) or 1:gt_L (females) for X-linked models
        #outputs vector of length |values| with genotype probs for indiv id given pedigree og partial genotype information
        probs = unlist(lapply(values, function(g) { 
            part = partial # to avoid NOTE in R CMD check
            part[id] = g  
            likelihood_LINKAGE(x, afreq=afreq, singleNum.geno=part, initialCalc=initialCalc, TR.MATR=.TRzero) 
        }))
        if (sum(probs)==0)  stop("\nIndividual ", x$orig.ids[id], ": All genotype probabilities zero. Mendelian error?") 
        probs
    }
    
    # do the rest of the individuals (only those present in sim_indivs)
    switch(chrom,    
    AUTOSOMAL = {
        for (i in cost_int) {
            markers[i, ] = apply(markers, 2, function(partgeno) sample.int(gt_L, size=1, prob=genoprobs(x, partial=partgeno, id=i, values=gt_values)) )
            likel_counter = likel_counter + length(gt_values)*N
        }
    },
    X = {
        for (i in setdiff(males_int, males_int[seq_len(init_m)])) {
            markers[i, ] = apply(markers, 2, function(partgeno) sample.int(gtM, size=1, prob=genoprobs(x, partial=partgeno, id=i, values=gt_values_m)) )
            likel_counter = likel_counter + length(gt_values_m)*N
        }
        for (i in setdiff(females_int, females_int[seq_len(init_f)])) {
            markers[i, ] = apply(markers, 2, function(partgeno) sample.int(gtF, size=1, prob=genoprobs(x, partial=partgeno, id=i, values=gt_values)) )
            likel_counter = likel_counter + length(gt_values)*N
        }
    })
   
    if(loops) { #quicker to tie loops before adding the markers
        markers = markers[-match(x$loop_breakers[, 2], x$orig.ids), , drop=F]
        x = tieLoops(x)
    }
   
    markers[!x$orig.ids %in% available, ] = 0
    if (unique) 
        markers = unique(markers, MARGIN=2)
    markers2 = .geno2diallelNEW(markers, nall)
    attrib = attributes(partialmarker)
    if(chrom == "X") attrib$chrom = 23
    markerdata_list = lapply(2*seq_len(ncol(markers)), function(k) {
        mk = markers2[, c(k-1, k)]
        attributes(mk) = attrib
        mk
    })
    class(markerdata_list) = "markerdata"
    x = setMarkers(x, markerdata_list)
    
    if(verbose) cat("\n", if(unique) x$nMark, if(unique) " unique markers simulated. ", likel_counter, " likelihood() calls.\nTotal time used: ", (proc.time() - starttime)[["elapsed"]], " seconds.\n", sep="")
    x
}
