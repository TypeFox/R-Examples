markerSim <- function(x, N=1, available=x$orig.ids, alleles=NULL, afreq=NULL, partialmarker=NULL, 
                     loop_breakers=NULL, eliminate=0, seed=NULL, method=3, verbose=TRUE) {
    starttime = proc.time()
    likel_counter = 0
    if (any(!is.numeric(N), length(N)>1, N%%1 != 0)) stop("N must be a positive integer.")
    if (!is.null(partialmarker)) {
        if (!inherits(partialmarker, "marker")) stop("Argument 'partialmarker' must be a 'marker' object.")
        else if (nrow(partialmarker)!=x$nInd) stop("Partial marker does not fit the pedigree.")
        if(!is.null(alleles) || !is.null(afreq)) stop("When 'partialmarker' is non-NULL, both 'alleles' and 'afreq' must be NULL.")
        if(length(mendelianCheck(setMarkers(x, partialmarker), verbose=F)) > 0) stop("Mendelian error in the given partial marker.")
    }
    else {
        if(is.null(alleles)) stop("Please specify marker alleles.")
        if(is.numeric(alleles) && length(alleles)==1) alleles = seq_len(alleles)
        partialmarker = marker(x, alleles=alleles, afreq=afreq)
    }
    m = partialmarker; alleles = attr(m, 'alleles'); afreq = attr(m, 'afreq')
    chrom = if(identical(23L, as.integer(attr(m, 'chrom')))) 'X' else 'AUTOSOMAL'

    if(all(m==0))  return(simpleSim(x, N, alleles=alleles, afreq=afreq, available=available, Xchrom=(chrom=='X'), seed=seed, verbose=verbose))

    allgenos = allGenotypes(nall <- attr(m, 'nalleles'))
    #al1 = allgenos[,1]; al2 = allgenos[,2]; nGeno = nrow(allgenos)
    
    if(verbose) {
        cat(ifelse(chrom=="AUTOSOMAL", "Autosomal", "X-linked"), "marker locus\n")
        cat("Simulating genotypes for", ifelse(length(available)==1, "individual", "individuals"), .prettycat(available, "and"), "\n")
        cat("\nAlleles and frequencies:\n")
        print(structure(round(afreq,3), names=alleles))
        cat("\nConditioning on the following genotypes:\n")
        print(data.frame(ID=x$orig.ids, GENO=.prettyMarkers(m, missing="-", singleCol=TRUE, sep="/", sex=x$pedigree[, 'SEX'])))
        cat("\n")
    }
    
    gridlist = geno.grid.subset(x, m, x$orig.ids, chrom, make.grid=F)
    for(id in 1:x$nInd)
        if(any(m[id, ]==0) && length(gridlist[[id]]) == 1) {
            m[id,] =  allgenos[gridlist[[id]], ]
            if(verbose) cat("Individual", x$orig.ids[id], "has forced genotype", ifelse(chrom=="X", alleles[m[id,1]], paste(alleles[m[id,]], collapse="/")), "\n")
        }

    preexisting = (m[,1]!=0 | m[,2]!=0)
    preexist_orig = x$orig.ids[preexisting]
    cost_orig = .mysetdiff(available, preexist_orig)
    if(method == 3) {  # TODO: Make this (a lot) smarter! Often worse than method 2, e.g. in twoloops pedigree. 
        sim_anc_orig = .mysetdiff(c(cost_orig, ancestors(x, id=cost_orig)), preexist_orig)
        cost_orig = sim_anc_orig[inpreanc <- (sim_anc_orig %in% ancestors(x, id=preexist_orig))]
        simpledrop_orig = sim_anc_orig[!inpreanc]
    } 
    
    if (loops <- x$hasLoops)    {        
        if(is.null(lb <- loop_breakers))      stop("The pedigree has loops. Please indicate loop breakers.")
        orig_ids = x$orig.ids
        if(verbose) cat(ifelse(length(lb)==1, "Breaking loop at individual ", "\nBreaking loops at individuals "), .prettycat(lb, "and"), "\n", sep="")
        x = breakLoops(setMarkers(x, partialmarker, missing=0), lb)
        m = x$markerdata[[1]]
        gridlist = gridlist[sort.int(match(c(orig_ids, lb), orig_ids))]
        if(method > 1)     cost_orig = unique.default(c(cost_orig, lb)) # added c() here.
    }

    cost_orig = cost_orig[order(!cost_orig %in% x$loop_breakers[,1])] #place loop breakers first!
    cost_int = .internalID(x, cost_orig)
    ngrid = lengths(gridlist)
    ped = x$pedigree    
    SEX = ped[,'SEX']
    # create initial marker matrix: two columns per marker
    markers = rep.int(m, N); dim(markers)=c(x$nInd, 2*N)
    odd = seq_len(N)*2 - 1
    if (!is.null(seed)) set.seed(seed)
    
    #Method = 2: pre-calculate probabilities for some individuals (big time saver!)
    if (method >= 2 && length(cost_orig)>0) {  
        init_int = switch(chrom,    
        AUTOSOMAL = {
            ngrid_cost = ngrid[cost_int]
            initvec = sapply(seq_along(cost_int), function(ci) 
                           prod(ngrid_cost[seq_len(ci)])  + N*sum(ngrid_cost[seq.int(ci+1, length.out=length(cost_int)-ci)]))
            cost_int[seq_len(which.min(initvec))]
        },
        X = {                
            males = cost_int[SEX[cost_int]==1]
            females = cost_int[SEX[cost_int]==2]
            ngrid_m = lengths(gridlist[males], use.names=F)
            ngrid_f = lengths(gridlist[females], use.names=F)
            n_males = length(males)
            n_females = length(females)
            # find optimal 'init' values for males/females
            calls = matrix(nrow=n_males+1, ncol=n_females+1)
            for(ma in 0:n_males) for(fe in 0:n_females)
                calls[ma+1, fe+1] = prod(ngrid_m[seq_len(ma)]) * prod(ngrid_f[seq_len(fe)]) + N*sum(c(ngrid_m[seq.int(ma+1, length.out=n_males-ma)], ngrid_f[seq.int(fe+1, length.out=n_females-fe)] )) # = number of times likelihood is called.

            calls.min = arrayInd(which.min(calls), dim(calls)) 
            c(males[seq_len(calls.min[1]-1)], females[seq_len(calls.min[2]-1)])
        })
        
        if(verbose) cat("\nTime saver: Pre-computing", ifelse(length(init_int)==1, "probabilities for individual", "joint probabilities for individuals"), .prettycat(sort(x$orig.ids[init_int]), 'and'), "\n")    
        
        allgenos_row_grid = t.default(fast.grid( gridlist[init_int] )) #Cartesian product. Each row contains 'init' row numbers of allgenos.
        initp = apply(allgenos_row_grid, 2, function(rownrs) { 
            partial = m
            partial[init_int, ] = allgenos[rownrs, ];   
            likelihood.linkdat(x, locus1=partial, eliminate=eliminate) 
        })
        likel_counter = likel_counter + length(initp)
        if (identical(sum(initp), 0)) stop("When trying to pre-compute joint probabilities: All probabilities zero. Mendelian error?")
        
        # fill the rows of the 'init' individuals 
        sample_rows = allgenos_row_grid[, suppressWarnings(sample.int(length(initp), size=N, replace=TRUE, prob=initp))]
        markers[init_int, odd] = allgenos[sample_rows, 1]
        markers[init_int, odd +1] = allgenos[sample_rows, 2]    
        cost_int = .mysetdiff(cost_int, init_int)
    }
    
    # The "costly" individuals. These require many likelihood calls: One call per marker per individual. (With 'method=1' everybody is treated this way.)
    if(verbose) {
        cat("\nBrute force (time consuming) conditional simulation ")
        if((lc <- length(cost_int))==0) cat("not needed\n") else cat("for", ifelse(lc==1, "individual", "individuals"), .prettycat(sort(x$orig.ids[cost_int]), 'and'), '\n')
    }
    for (i in cost_int) {
        gridi = gridlist[[i]]
        rowsample = unlist(lapply(2*seq_len(N), function(mi) {
            partial = m
            partial[] = markers[, c(mi-1, mi)]  # preserves all attributes of the m.
            probs = unlist(lapply(gridi, function(r) { partial[i, ] = allgenos[r,];    likelihood.linkdat(x, locus1=partial, eliminate=eliminate)    }))
            if (sum(probs)==0) { print(cbind(ped, partial)); stop("\nIndividual ", x$orig.ids[i],": All genotype probabilities zero. Mendelian error?")}
            sample(gridi, size=1, prob=probs)
        }))
        markers[i, odd] = allgenos[rowsample, 1]; markers[i, odd+1] = allgenos[rowsample, 2]  
    }
    likel_counter = likel_counter + N*sum(ngrid[cost_int])
    
    # Method=3: Simulate final individuals by sampling random founder alleles followed by gene dropping:
    if(method == 3 && length(simpledrop_orig)>0) {
        if(verbose) cat("\nTime saver: Simulation by gene dropping for individuals", .prettycat(sort(simpledrop_orig), 'and'), '\n')
        
        loopbr_int = .internalID(x, x$loop_breakers[, 1]) #integer(0) if no loops
        loopbr_dup_int = .internalID(x, x$loop_breakers[, 2])
        
        simpledrop_int = .internalID(x, simpledrop_orig)
        simple_founders = simpledrop_int[simple_is_founder <- simpledrop_int %in% x$founders]
        simple_nonfounders = sort.int(simpledrop_int[!simple_is_founder]) #sorting is crucial here, to make sure genedropping goes down the pedigree.
        
        if(chrom=='AUTOSOMAL')     markers[simple_founders, ] = sample.int(nall, size=2*N*sum(simple_is_founder), replace=TRUE, prob=afreq)
        else for (f in simple_founders)  
            markers[f, ] = switch(SEX[f], rep(sample.int(nall, size=N, replace=TRUE, prob=afreq), each=2), sample.int(nall, size=2*N, replace=TRUE, prob=afreq))
        
        markers[loopbr_dup_int,] = markers[loopbr_int, ] # Genotypes of the duplicated individuals. Some of these may be ungenotyped...save time by excluding these?
        
        for(id in simple_nonfounders) {
            fa = ped[id, 'FID']; mo = ped[id, 'MID']
            if(chrom=='AUTOSOMAL') {
                markers[id, odd] = markers[fa, odd + .rand01(N)];    markers[id, odd+1] = markers[mo, odd + .rand01(N)]
            }
            else {
                switch(SEX[id], { markers[id, ] <- rep(markers[mo, odd + .rand01(N)], each=2) },
                    {markers[id, odd] <- markers[fa, odd + .rand01(N)];        markers[id, odd+1] <- markers[mo, odd + .rand01(N)] }
            )}    
            if((indx <- match(id, loopbr_int, nomatch=0)) > 0)
                markers[loopbr_dup_int[indx], ] = markers[id, ]                
        }
    }
            
    if(loops) {
        markers = markers[-match(x$loop_breakers[,2], x$orig.ids), ]
        x = tieLoops(x)
    }
    markers[!preexisting & !(x$orig.ids %in% available), ] = 0
    attrib = attributes(partialmarker)
    markerdata_list = lapply(seq_len(N), function(k) {mk = markers[, c(2*k-1,2*k)];    attributes(mk) = attrib; mk})
    class(markerdata_list) = "markerdata"
    x = setMarkers(x, markerdata_list)
    if(verbose) cat("\n",x$nMark, " markers simulated.\nNumber of calls to the likelihood function: ", likel_counter, ".\nTotal time used: ", (proc.time() - starttime)[["elapsed"]], " seconds.\n", sep="")
    x
}


simpleSim = function(x, N, alleles, afreq, available, Xchrom=FALSE, seed=NULL, verbose=T) {
    starttime = proc.time()
    if(missing(alleles)) {
        if(missing(afreq)) stop("Both 'alleles' and 'afreq' cannot be missing")
        alleles = seq_along(afreq)
    }
    nall = length(alleles)
    if(missing(afreq)) 
        afreq = rep(1, nall)/nall
    if(variableSNPfreqs <- (nall==2 && length(afreq)!=2 && !Xchrom))
        afreq = rep(afreq, length=N)
    if(missing(available)) 
        available = x$orig.ids
        
    if(verbose) {
        cat(ifelse(!Xchrom, "Autosomal", "X-linked"), "marker locus\n")
        cat("Unconditional simulating of genotypes for", ifelse(length(available)==1, "individual", "individuals"), .prettycat(available, "and"), "\n")
        cat("\nAlleles and frequencies:\n")
        if(variableSNPfreqs) cat("SNPs with allele 1 frequencies", paste(head(afreq, 5), collapse=", "), ifelse(N>5, '...\n','\n'))
        else print(structure(round(afreq,3), names=alleles))
    }
    
    ped = x$pedigree
    m = matrix(0, ncol=2*N, nrow=x$nInd)
    odd = seq_len(N)*2 - 1
    
    if (!is.null(seed)) set.seed(seed)
    if(Xchrom) {
        for (f in x$founders) {
            if(ped[f, 'SEX']==1) m[f, ] = rep(sample.int(nall, size=N, replace=TRUE, prob=afreq), each=2)
            if(ped[f, 'SEX']==2) m[f, ] = sample.int(nall, size=2*N, replace=TRUE, prob=afreq)
        }
        for(id in x$nonfounders) {
            fa = ped[id, 'FID']; mo = ped[id, 'MID']
            switch(ped[id, 'SEX'], 
                {m[id, ] <- rep(m[mo, odd + .rand01(N)], each=2)},
                {m[id, odd] <- m[fa, odd + .rand01(N)]; m[id, odd+1] = m[mo, odd + .rand01(N)]})
        }
    }
    else {
        size = 2*length(x$founders)
        allelsamp = if(variableSNPfreqs) unlist(lapply(afreq, function(f) sample.int(2, size, replace=TRUE, prob=c(f,1-f))))
                    else sample.int(nall, size=N*size, replace=TRUE, prob=afreq)
        m[x$founders, ] = allelsamp
        for(id in x$nonfounders) {
            fa = ped[id, 'FID']; mo = ped[id, 'MID']
            m[id, odd] = m[fa, odd + .rand01(N)]
            m[id, odd+1] = m[mo, odd + .rand01(N)]
        }
    }

    m[!x$orig.ids %in% available, ] = 0
    if(variableSNPfreqs) {
        attrib = attributes(marker(x, alleles=alleles, afreq=NULL, chrom=NA, missing=0))
        frqs = as.vector(rbind(afreq,1-afreq))
        markerdata_list = lapply(odd, function(k) {
            mk = m[, c(k, k+1)]
            atr = attrib; atr$afreq = frqs[c(k, k+1)]; 
            attributes(mk) = atr
            mk
        })
    }
    else {
        attrib = attributes(marker(x, alleles=alleles, afreq=afreq, chrom=ifelse(Xchrom, 23, NA), missing=0))
        markerdata_list = lapply(odd, function(k) {mk = m[, c(k, k+1)];    attributes(mk) = attrib; mk})
    }
    x = setMarkers(x, structure(markerdata_list, class = "markerdata"))
    if(verbose) cat("\n", x$nMark, " markers simulated.\nNumber of calls to the likelihood function: 0.\nTotal time used: ", (proc.time() - starttime)[["elapsed"]], " seconds.\n", sep="")
    x
}


