twoLocusIBD = function(x, ind1, ind2, rho=NULL, cM=NULL, Nsim, Xchrom=FALSE, verbose=TRUE, ...) {
    if((err<-ind1) > x$nInd || (err<-ind2) > x$nInd) stop(paste("The pedigree has no individual with label", err))
    if(is.null(cM) + is.null(rho) != 1) stop("Exactly one of the parameters 'cM' and 'rho' must be non-NULL.")
    starttime = proc.time()
    
    if(is.null(cM)) {
        if(rho<0 | rho>0.5) stop("Recombination rate 'rho' is outside of interval [0, 0.5].")
        cM = -50*log(1-2*rho)
    }
    
    if(cM==Inf) {
        if(verbose) cat("Analysis of unlinked loci.\n")
        m1 = oneLocusIBD(x, ind1, ind2, Nsim=Nsim, Xchrom=Xchrom, verbose=verbose, ...)
        m2 = oneLocusIBD(x, ind1, ind2, Nsim=Nsim, Xchrom=Xchrom, verbose=verbose)
        res = outer(m1, m2)
        return(res)
    }
    if(verbose) cat("Locus distance:", cM, "centiMorgan\n")
    map = uniformMap(cM = cM, chromosome=if(Xchrom) 23 else 1)
    simdata = IBDsim(x, map=map, sims=Nsim, model="haldane", verbose=verbose, ...)
    
    # Utility function: IBD status for a pair of (non-inbred) genotypes. Each genotype is a pair of alleles.
    ibd.status = function(gt1, gt2)
        (gt1[1] %in% gt2) + (gt1[2] %in% gt2)
        
    # Setup for X chromosomal analysis
    if(Xchrom) {
        sex1 = x$pedigree[ind1, "SEX"] # NB: IBDsim raises error if not labels are 1,2,3,...
        sex2 = x$pedigree[ind2, "SEX"]
        
        if(sex1==2 && sex2==1) { # reverse indiv order in this case (simplifies code below)
            inds = c(ind1, ind2)
            ind1 = inds[2]; ind2=inds[1]
            sex1 = 1; sex2 = 2
        }
    }
    
    # Define appropriate function for analysing a single simulation: Output is a pair (IBD_locus1, IBD_locus2)
    if(Xchrom && sex1==1 && sex2==1) 
        single.sim.ibd = function(sim) {
            i1.mat = sim[[c(1, ind1, 2)]] # maternal chromosome of indiv 1 (matrix with 2 columns: breakpoint position; allele)
            i2.mat = sim[[c(1, ind2, 2)]]
            
            # marker 1
            m1_ibd = i1.mat[1, 2] == i2.mat[1, 2] # first row (first locus), second column (allele)
            
            # marker 2
            m2_ibd =  i1.mat[dim(i1.mat)[1], 2] == i2.mat[dim(i2.mat)[1], 2] #last row (second locus), second column (allele)
            c(m1_ibd, m2_ibd)
        }
    else if(Xchrom && sex1==1) # includes sex1=2, sex2=1, because of swapping above
        single.sim.ibd = function(sim) {
            i1.mat = sim[[c(1, ind1, 2)]] # maternal chromosome of indiv 1 (matrix with 2 columns: breakpoint position; allele)
            i2.pat = sim[[c(1, ind2, 1)]] 
            i2.mat = sim[[c(1, ind2, 2)]]
            
            # marker 1
            ind1.mat.m1 = i1.mat[1, 2]
            ind2.pat.m1 = i2.pat[1, 2]
            ind2.mat.m1 = i2.mat[1, 2]
            m1_ibd =  ind1.mat.m1 %in% c(ind2.pat.m1, ind2.mat.m1)
            
            # marker 2
            ind1.mat.m2 = i1.mat[dim(i1.mat)[1], 2] #last row (second locus), second column (allele)
            ind2.pat.m2 = i2.pat[dim(i2.pat)[1], 2]
            ind2.mat.m2 = i2.mat[dim(i2.mat)[1], 2]
            m2_ibd =  ind1.mat.m2 %in% c(ind2.pat.m2, ind2.mat.m2)
            
            c(m1_ibd, m2_ibd)
        }
    else 
        single.sim.ibd = function(sim) {
            i1.pat = sim[[c(1, ind1, 1)]] # paternal chromosome of indiv 1 (matrix with 2 columns: breakpoint position; allele)
            i1.mat = sim[[c(1, ind1, 2)]]
            i2.pat = sim[[c(1, ind2, 1)]] 
            i2.mat = sim[[c(1, ind2, 2)]]
            
            # marker 1
            ind1.pat.m1 = i1.pat[1, 2]
            ind1.mat.m1 = i1.mat[1, 2]
            ind2.pat.m1 = i2.pat[1, 2]
            ind2.mat.m1 = i2.mat[1, 2]
            m1_ibd =  ibd.status(c(ind1.pat.m1, ind1.mat.m1), c(ind2.pat.m1, ind2.mat.m1))
            
            # marker 2
            ind1.pat.m2 = i1.pat[dim(i1.pat)[1], 2]
            ind1.mat.m2 = i1.mat[dim(i1.mat)[1], 2]
            ind2.pat.m2 = i2.pat[dim(i2.pat)[1], 2]
            ind2.mat.m2 = i2.mat[dim(i2.mat)[1], 2]
            m2_ibd =  ibd.status(c(ind1.pat.m2, ind1.mat.m2), c(ind2.pat.m2, ind2.mat.m2))
            
            c(m1_ibd, m2_ibd)
        }
    
    
    # Analyse simulations
    allsims.ibd = vapply(simdata, single.sim.ibd, numeric(2))
    
    # Frequency table of results
    max.ibd = if(Xchrom && sex1==1) 1 else 2
    res = matrix(ncol = max.ibd+1, nrow = max.ibd+1, dimnames = rep(list(paste0("ibd", 0:max.ibd)), 2))
    for(i in 0:max.ibd) for(j in 0:max.ibd)
        res[i+1,j+1] = sum((allsims.ibd[1,] == i) & (allsims.ibd[2,] == j))
    
    if(verbose) cat("IBD analysis finished. Total time used:", (proc.time()-starttime)[["elapsed"]], "seconds.\n")
    res/Nsim
}


oneLocusIBD = function(x, ind1, ind2, Nsim, Xchrom=FALSE, verbose=TRUE, ...) {
    if((err<-ind1) > x$nInd || (err<-ind2) > x$nInd) stop(paste("The pedigree has no individual with label", err))
    
    # Utility function: IBD status for a pair of (non-inbred) genotypes. Each genotype is a pair of alleles.
    ibd.status = function(gt1, gt2)
        (gt1[1] %in% gt2) + (gt1[2] %in% gt2)
        
    # Setup for X chromosomal analysis
    if(Xchrom) {
        sex1 = x$pedigree[ind1, "SEX"] # NB: IBDsim raises error if not labels are 1,2,3,...
        sex2 = x$pedigree[ind2, "SEX"]
        
        if(sex1==2 && sex2==1) { # reverse indiv order in this case (simplifies code below)
            inds = c(ind1, ind2)
            ind1 = inds[2]; ind2=inds[1]
            sex1 = 1; sex2 = 2
        }
    }
    
    # Define appropriate function for analysing a single simulation: Output is a pair (IBD_locus1, IBD_locus2)
    if(Xchrom && sex1==1 && sex2==1) 
        single.sim.ibd = function(sim) {
            i1.mat = sim[[c(1, ind1, 2)]] # maternal chromosome of indiv 1. 1*2 matrix: [0, allele]
            i2.mat = sim[[c(1, ind2, 2)]]
            i1.mat[1, 2] == i2.mat[1, 2] # TRUE iff IBD = 1
        }
    else if(Xchrom && sex1==1) # includes sex1=2, sex2=1, because of swapping above
        single.sim.ibd = function(sim) {
            i1.mat = sim[[c(1, ind1, 2)]] # maternal chromosome of indiv 1
            i2.pat = sim[[c(1, ind2, 1)]] 
            i2.mat = sim[[c(1, ind2, 2)]]
            # each of the above is a 1*2 matrix: [0, allele]
            
            i1.mat[1, 2] %in% c(i2.pat[1, 2], i2.mat[1, 2])
        }
    else 
        single.sim.ibd = function(sim) {
            i1.pat = sim[[c(1, ind1, 1)]] # paternal chromosome of indiv 1 (matrix with 2 columns: breakpoint position; allele)
            i1.mat = sim[[c(1, ind1, 2)]]
            i2.pat = sim[[c(1, ind2, 1)]] 
            i2.mat = sim[[c(1, ind2, 2)]]
            # each of the above is a 1*2 matrix: [0, allele]
            
            ibd.status(c(i1.pat[1, 2], i1.mat[1, 2]), c(i2.pat[1, 2], i2.mat[1, 2]))
        }
    
    starttime = proc.time()
    
    map = uniformMap(cM = 0, chromosome=if(Xchrom) 23 else 1)
    simdata = IBDsim(x, map=map, sims=Nsim, model="haldane", verbose=verbose, ...)
    
    # Analyse simulations
    allsims.ibd = vapply(simdata, single.sim.ibd, numeric(1))
    
    ibdlevels = if(Xchrom && sex1==1) 0:1 else 0:2
    res = table(factor(allsims.ibd, levels=ibdlevels, labels=paste0('ibd', ibdlevels)))
    res = c(res) # strips table attributes, but keeps names. (Hadley: Bad style.)
    if(verbose) cat("IBD analysis finished. Total time used:", (proc.time()-starttime)[["elapsed"]], "seconds.\n")
    res/Nsim
}

