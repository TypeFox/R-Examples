# expectedGenotypes.R
#--------------------------------------------------------------------------
# What: Construct expected genotypes according to known allele variants
# Time-stamp: <2007-07-19 12:06:26 ggorjan>
#--------------------------------------------------------------------------

expectedGenotypes <- function(x, alleles=allele.names(x), ploidy=2,
                              sort=TRUE, haplotype=FALSE)
{
    ## Checks
    if (missing(x) && missing(alleles))
        stop("at least one of 'x' or 'alleles' must be given")
    if (!(missing(x) && !missing(alleles))) {
        if (!is.genotype(x)) stop("x must be of class 'genotype'")
    }

    nA <- length(alleles)
    if(nA==0)
      stop("Zero length alleles vector.")

    ## Add numbers to allele names, so we get result sorted by
    ## given allele names
    if (sort) {
        allelesOrig <- alleles
        alleles <- paste(1:nA, alleles, sep="")
    }

    ## Find possible genotypes according to allele variants
    if(!haplotype) {
      comb <- combinations(n=nA, r=ploidy, v=alleles, repeats.allowed=T)
    } else {
      comb <- permutations(n=nA, r=ploidy, v=alleles, repeats.allowed=T)
    }

    if (sort) {
        for (i in 1:nA) comb[comb == alleles[i]] <- allelesOrig[i]
    }

    ## Create a nice character vector of expected genotypes
    nC <- dim(comb)
    ret <- vector(mode="character", length=nC[1])
    ret <- paste(comb[, 1], comb[, 2], sep="/")
    if (ploidy > 2) {
        for (i in 3:nC[2]) {
            ret <- paste(ret, comb[, i], sep="/")
        }
    }
    return(ret)
}

expectedHaplotypes <- function(x, alleles=allele.names(x), ploidy=2,
                              sort=TRUE, haplotype=TRUE)
{
  expectedGenotypes(x=x, alleles=alleles, ploidy=ploidy,
                    sort=sort, haplotype=haplotype)
}

#--------------------------------------------------------------------------
# expectedGenotypes.R ends here
