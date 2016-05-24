
twoMarkerDistribution <- function(x, id, partialmarker1, partialmarker2, theta, loop_breakers=NULL, eliminate=99, verbose=TRUE) {
    starttime = proc.time()
    if (!inherits(m1 <- partialmarker1, "marker"))
        if(is.numeric(m1) && length(m1)==1 && m1 <= x$nMark) m1 = x$markerdata[[m1]] 
        else stop("The 'partialmarker1' must be a 'marker' object, or a single integer indicating an existing marker of 'x'.")
    if (!inherits(m2 <- partialmarker2, "marker"))
        if(is.numeric(m2) && length(m2)==1 && m2 <= x$nMark) m2 = x$markerdata[[m2]] 
        else stop("The 'partialmarker2' must be a 'marker' object, or a single integer indicating an existing marker of 'x'.")
    
    markerchrom1 = as.integer(attr(m1, 'chrom'))
    markerchrom2 = as.integer(attr(m2, 'chrom'))
    if(!identical(markerchrom1, markerchrom2)) 
        stop(sprintf("Partial markers are on different chromosomes: %d and %d.", markerchrom1, markerchrom2))
    chrom = ifelse(identical(23L, markerchrom1), 'X', 'AUTOSOMAL')
    
    SEX = x$pedigree[,'SEX']   
    afreq1 = attr(m1, 'afreq'); alleles1 = attr(m1, 'alleles')
    afreq2 = attr(m2, 'afreq'); alleles2 = attr(m2, 'alleles')
    
    if(verbose) {
        cat(ifelse(chrom=="AUTOSOMAL", "Autosomal", "X-linked"), "markers with the following known genotypes:\n")
        print(data.frame(ID=x$orig.ids, 
                         M1=as.character(.prettyMarkers(list(m1), missing="-", singleCol=TRUE, sex=SEX)),
                         M2=as.character(.prettyMarkers(list(m2), missing="-", singleCol=TRUE, sex=SEX))), row.names=FALSE)
        cat("\nAllele frequencies, marker 1:\n"); print(structure(afreq1, names=alleles1))
        cat("\nAllele frequencies, marker 2:\n"); print(structure(afreq2, names=alleles2))
        cat("\nRecombination rate between marker loci:", theta,"\n")
    }
    
    #Do this before loop breaking, since eliminate2 works better WITH the loops.
    grid.subset = fast.grid(c(geno.grid.subset(x, m1, id, chrom, make.grid=F), 
                              geno.grid.subset(x, m2, id, chrom, make.grid=F)))
    
    if (x$hasLoops) {
        if(is.null(lb <- loop_breakers))      stop("The pedigree has loops. Please indicate loop breakers.")
        if(verbose) cat(ifelse(length(lb)==1, "\nBreaking loop at individual", "\nBreaking loops at individuals"), .prettycat(lb, "and"), "\n")
        x = breakLoops(setMarkers(x, list(m1, m2)), lb)
        m1 = x$markerdata[[1]]
        m2 = x$markerdata[[2]]
        SEX = x$pedigree[,'SEX']
    }
    
    int.id = .internalID(x, id)
    allgenos1 = allGenotypes(attr(m1, 'nalleles'))
    allgenos2 = allGenotypes(attr(m2, 'nalleles'))
    
    if(chrom=='AUTOSOMAL' || SEX[int.id]==2) {
        gt1.strings = paste(alleles1[allgenos1[, 1]], alleles1[allgenos1[, 2]], sep="/")
        gt2.strings = paste(alleles2[allgenos2[, 1]], alleles2[allgenos2[, 2]], sep="/")
        geno.names = list(gt1.strings, gt2.strings)
    } else
        geno.names = list(alleles1, alleles2)
    
    marginal = likelihood(x, locus1=m1, locus2=m2, theta=theta, eliminate=eliminate)
    if(marginal==0) stop("Partial marker data is impossible")
    
    probs = array(0, dim=lengths(geno.names, use.names=F), dimnames=geno.names)
    probs[grid.subset] = apply(grid.subset, 1, function(allg_rows) {
        m1[int.id, ] = allgenos1[allg_rows[1], ]
        m2[int.id, ] = allgenos2[allg_rows[2], ]
        likelihood(x, locus1=m1, locus2=m2, theta=theta, eliminate=eliminate)
    })
    
    res = probs/marginal
    if(verbose) {
        cat("\nJoint genotype distribution at the two markers for individual ", id,":\n", sep="")
        print(round(res,4))
        cat("\nTotal time used: ", (proc.time() - starttime)[["elapsed"]], " seconds.\n", sep="")
        return(invisible(res))
    }
    else res
}
