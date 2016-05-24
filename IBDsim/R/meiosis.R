meiosis <-
function(parent, map, model="chi", condition=NULL, skip.recomb=FALSE) { #skip=TRUE returns random strand with no recombination; condition should be NULL or a vector with elements 'locus'(Mb), 'allele' and 'action' (1=force,2=avoid).
    if (condit <- !is.null(condition)) {
        whichStrand = which(condition[['allele']] == .getAlleles(parent, locus <- condition[['locus']]))
        startStrand = switch(condition[['action']], 
            switch(length(whichStrand) + 1, stop("Allele to be forced is not present."), whichStrand, sample.int(2,1)),
            switch(length(whichStrand) + 1, sample.int(2,1), 3 - whichStrand, stop("Allele cannot be avoided."))
        )
    } else     startStrand = sample.int(2,1)
    
    if (skip.recomb) return(parent[[startStrand]])
    L.cM = map[nrow(map), 'cM']  #chromosome length in cM
    
    switch(model, 
    haldane = {
        ncross = as.integer(rpois(1, L.cM/100))
        if (ncross==0) return(parent[[startStrand]])
        Cx = .sortDouble(runif(ncross, min=0, max=L.cM))
    },
    chi = {
        m = 4
        nC = rpois(1, L.cM/50*(m+1))    #L.cM/100*2*(m+1); number of potential crossover events
        if (nC==0) return(parent[[startStrand]])
        C_events = .sortDouble(runif(nC, min=0, max=L.cM)) #potential crossover positions (N-1 intervals, uniformly distr given nC)
        Cx.bundle = C_events[!as.logical((seq_len(nC) + sample.int(m + 1, 1)) %% (m + 1))]    #Cx events on 4 strand bundle: every (m+1)th
        Cx = Cx.bundle[as.logical(sample.int(2, length(Cx.bundle), replace=T)%%2)]    #thinning. Each survive with prob=1/2
        ncross = length(Cx)
        if (ncross==0) return(parent[[startStrand]])
    })
    cpos = cm2phys(cM_locus = Cx, mapmat=map) #crossover positions
    if(condit)     startStrand = 2 - (startStrand + sum(cpos < locus))%%2  #switches start strand iff sum(cpos < locus) is odd
        
    par1 = parent[[startStrand]]; p1 = as.double(par1[,1]); a1 = as.integer(par1[,2]); l1=length(p1)
    par2 = parent[[3-startStrand]]; p2 = as.double(par2[,1]); a2 = as.integer(par2[,2]); l2=length(p2)
    pos = numeric(l1+l2+ncross)
    allel = integer(l1+l2+ncross)
    resC = .C("recombine", p1, a1, l1, p2, a2, l2, cpos, ncross, pos=pos, allel=allel, PACKAGE="IBDsim")
    keep = seq_len(match(0, resC$pos[-1]))
    res = c(resC$pos[keep], resC$allel[keep])
    dim(res) = c(length(res)/2,2)
    res    
}


.sortDouble = function(x) x[order(x)] #.Internal(qsort(x, FALSE))
