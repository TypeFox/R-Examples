
oneMarkerDistribution <- function(x, ids, partialmarker, theta=NULL, grid.subset=NULL, loop_breakers=NULL, eliminate=0, verbose=TRUE) {
    starttime = proc.time()
    if (!inherits(m <- partialmarker, "marker"))
        if(is.numeric(m) && length(m)==1 && m <= x$nMark) 
            m = x$markerdata[[m]] 
        else stop("The 'partialmarker' must be a 'marker' object, or a single integer indicating an existing marker of 'x'.")
    markerchrom = as.integer(attr(m, 'chrom'))
    alleles = attr(m, "alleles")
    
    if (affped <- any(x$pedigree[, 'AFF'] == 2)) {
        if(is.null(theta)) stop("The pedigree is affected with disease: Please spesify the recombination fraction ('theta') between the marker and disease loci.")
        if(is.null(x$model)) stop("The pedigree is affected with disease: Please spesify a disease model, e.g. using 'setModel()'.")
        locus2 = 'disease'
        chrom = x$model$chrom
        if(!is.na(markerchrom) && ((chrom == 'X') != (23L == markerchrom))) stop("Disease model and marker chromosome are not compatible.")
    } 
    else {
        locus2 = NULL
        chrom = ifelse(identical(23L, markerchrom), 'X', 'AUTOSOMAL')
    }
        
    SEX = x$pedigree[,'SEX']
    if(verbose) {
        cat(ifelse(chrom=="AUTOSOMAL", "Autosomal", "X-linked"), "marker with the following partial data:\n")
        print(data.frame(ID=x$orig.ids, GENO=.prettyMarkers(list(m), missing="-", singleCol=TRUE, sex=SEX)), row.names=FALSE)
        cat("\nMarker allele frequencies:\n")
        print(structure(attr(m, 'afreq'), names=alleles))
        if(affped) {cat("\nDisease model:\n"); print(x$model);    cat("\nRecombination rate between marker and disease locus: ", theta,".\n", sep="")}
    }
    
    allgenos = allGenotypes(attr(m, 'nalleles'))
    
    if(is.null(grid.subset)) grid.subset = geno.grid.subset(x, m, ids, chrom) #Do this before loop breaking, works better with eliminate2.
    else grid.subset = as.matrix(grid.subset) 
    
    if (x$hasLoops) {
        if(is.null(lb <- loop_breakers))      stop("The pedigree has loops. Please indicate loop breakers.")
        if(verbose) cat(ifelse(length(lb)==1, "\nBreaking loop at individual", "\nBreaking loops at individuals"), .prettycat(lb, "and"), "\n")
        x = breakLoops(setMarkers(x, m), lb)
        m = x$markerdata[[1]]
        SEX = x$pedigree[,'SEX']
    }
    
    int.ids = .internalID(x, ids)
    gt.strings = paste(alleles[allgenos[, 1]], alleles[allgenos[, 2]], sep="/")
    geno.names = switch(chrom, 
        AUTOSOMAL = rep(list(gt.strings), length(ids)),
        X = list(alleles, gt.strings)[SEX[int.ids]]
    )
    
    marginal = likelihood(x, locus1=m, locus2=locus2, theta=theta, eliminate=eliminate)
    if(marginal==0) stop("Partial marker is impossible")
    probs = array(0, dim=lengths(geno.names, use.names=F), dimnames=geno.names)
    probs[grid.subset] = apply(grid.subset, 1, function(allg_rows) {
        m[int.ids, ] = allgenos[allg_rows, ]
        likelihood(x, locus1=m, locus2=locus2, theta=theta, eliminate=eliminate)
    })
    
    res = probs/marginal
    if(verbose) {
        cat(ifelse(length(ids)==1, "\nGenotype probability distribution for individual ",
                "\nJoint genotype probability distribution for individuals "), .prettycat(ids,"and"), ":\n", sep="")
        print(round(res,4))
        cat("\nTotal time used: ", (proc.time() - starttime)[["elapsed"]], " seconds.\n", sep="")
        return(invisible(res))
    }
    else res
}
