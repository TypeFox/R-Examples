IBDsim <-
function(x, sims, query=NULL, condition=NULL, map="decode", chromosomes=NULL, 
    model="chi", merged=TRUE, simdata=NULL, skip.recomb = "noninf_founders", seed=NULL, verbose=TRUE) {
    
    if(!all(x$orig.ids==1:x$nInd)) stop("Individual ID's must be 1, 2, 3, ... . Please relabel (see e.g. ?relabel).")
    starttime = proc.time()
        
    if (is.null(simdata)) {
        if(!is.null(seed)) set.seed(seed)
        map = loadMap(map, chrom = chromosomes)
        if (verbose) {
            cat("---------------\n")
            cond = !is.null(condition)
            cat("Performing", ifelse(cond, 'conditional', 'unconditional'), "simulation\n") 
            mapchrom = attr(map, 'chromosome')
            if(is.null(mapchrom)) mapchrom = sapply(map, attr, 'chromosome')
            cat('Mode:', ifelse(identical(mapchrom, 23), "X chromosome", ifelse(23 %in% mapchrom, 'both', 'autosomal')), '\n')
            cat('Recombination model:', match.arg(model, c("haldane (poisson process)", "chi square renewal model")), '\n')
            cat('Number of simulations:', sims, '\n')
        }   
        
        if(!is.null(condition)) {
            if(length(map)==1) dischr = rep.int(attr(map[[1]], 'chromosome'), sims)
            else dischr = sample(sapply(map, attr, 'chromosome'), size=sims, replace=T, prob=sapply(map, attr, 'length_Mb'))
            oblig.saps = sample.obligates(x, condition, sims)
        }
        else {dischr=rep.int(0, sims)}
    
        if(!is.null(skip.recomb)) {
            if(skip.recomb=="noninf_founders") {
                cafs = x$founders; if(!is.null(condition)) cafs = intersect(cafs, .CAFs(x,condition)); if(!is.null(query)) cafs = intersect(cafs, .CAFs(x,query))
                skip.recomb = setdiff(x$founders, cafs)
            }
            if (length(skip.recomb)>0 && verbose) cat("Skipping recombination in the following individuals:", paste(skip.recomb, collapse=", "),"\n")
        }

        simdata = lapply(1:sims, function(i) 
            lapply(map, function(m) {
                if(dischr[sims]==attr(m, 'chromosome')) cond=oblig.saps[[i]] else cond=NULL
                genedrop(x, map=m, condition=cond, model=model, skip.recomb=skip.recomb)
            }))
        attr(simdata, 'total_map_length_Mb') = attr(map, "length_Mb")
        if (verbose) cat("Simulation finished. Time used:", 
                         (proc.time()-starttime)[['elapsed']], "seconds\n---------------\n")
    }
    if(is.null(query)) return(invisible(simdata))
    
    coeffs = inbreeding(x); inbreds = which(coeffs>0); 
    if(length(inbreds)>0) {
        inb = cbind(ID=inbreds, f=coeffs[inbreds]); rownames(inb) = rep("", length(inbreds))
        if (verbose) {cat("\nInbreeding coefficients:\n"); print(inb)}
    }
        
    runs <- lapply(simdata, function(h)     sap.segments(h, sap=query))
    attr(runs, 'total_map_length_Mb') = attr(simdata, 'total_map_length_Mb')
    
    if (verbose) cat("\nResults:\n")
    stats = summary.ibd(runs, merged=merged, verbose=verbose)
    if (verbose) cat("\nTotal time used:", (proc.time()-starttime)[['elapsed']], "seconds.\n")
    
    invisible(list(simdata=simdata, segments=runs, stats=stats))
}

sample.obligates = function(x, condition, sims) {
    obligate_ones = obligate.carriers(x, condition)
    complete.saps = lapply(obligate_ones, function(ones) {sap = condition; sap[['1']] = ones; sap})
    if(length(complete.saps)==1) {cat("For the disease chromosome I'm conditioning on the following SAP:\n"); .printSAP(complete.saps[[1]])}
    else {
        cat("For the disease chromosome I'm sampling condition SAPs among the following:\n") 
        for(i in 1:length(complete.saps)) {cat("SAP ",i,":\n",sep=""); .printSAP(complete.saps[[i]])}
    }
    weight = sapply(obligate_ones, function(vec) .5^(length(vec)-1))
    oblig.samples = sample(complete.saps, size=sims, replace=TRUE, prob=weight)
}


.printSAP = function(sap) {
    if(!is.null(two <- sap[['2']])) cat("  Two copies:", paste(two, collapse=", "), "\n")
    if(!is.null(atl1 <- sap[['atleast1']])) cat("  At least one copy:", paste(atl1, collapse=", "), "\n")
    if(!is.null(one <- sap[['1']])) cat("  One copy:", paste(one, collapse=", "), "\n")
    if(!is.null(atm1 <- sap[['atmost1']])) cat("  At most one copy:", paste(atm1, collapse=", "), "\n")
    if(!is.null(zero <- sap[['0']])) cat("  Zero copies:", paste(zero, collapse=", "), "\n")
    cat("\n")
}

.prettycat = function(v, andor)
    switch(min(len <- length(v), 3), toString(v), paste(v, collapse=" and "), paste(paste(v[-len], collapse=", "), andor, v[len]))