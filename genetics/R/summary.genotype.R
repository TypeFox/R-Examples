# $Id: summary.genotype.R 1314 2007-09-12 10:41:11Z ggorjan $

###
### Provide the frequency and proportions of alleles and genotypes
###

# used when summary.genotype is called from summary.data.frame:
shortsummary.genotype <- function(object, ..., maxsum)
  {
    tmp <- summary.factor(object, maxsum=maxsum)
    retval <- paste(format(tmp), " (", format(round(prop.table(tmp)*100)), "%)", sep='' )
    names(retval) <- names(tmp)
    #retval <- retval[order(tmp, decreasing=TRUE)]
    retval
  }

# general function
summary.genotype  <-  function(object,...,maxsum)
  {
    # if we are called from within summary.data.frame, fall back to
    # summary.factor so that we don't mess up the display
    if(!missing(maxsum))
      return(shortsummary.genotype(object,...,maxsum=maxsum))

    retval  <-  list()
#    retval$genotype  <- object
    retval$allele.names  <- allele.names(object)

    retval$locus  <- attr(object,"locus")
    class(retval)  <- "summary.genotype"

    af  <- table(allele(object))

    # make sure af has same order as allele.names...
    #
    missed <- !names(af) %in% retval$allele.names 
    af.tab <- rep(0,length(retval$allele.names))
    names(af.tab) <- retval$allele.names
    af.tab[names(af)] <- af
    #

    paf <- prop.table(af.tab)
    retval$allele.freq    <- cbind("Count"=af.tab,"Proportion"=paf)

    gf  <- table( object )
    pgf <- prop.table(gf)
    retval$genotype.freq  <- cbind("Count"=gf,"Proportion"=pgf)

    ## Sort by genotypeOrder
    asFun <- as.genotype
    if(is.haplotype(object)) asFun <- as.haplotype
    tmp <- asFun(rownames(retval$genotype.freq),
                 alleles=allele.names(object))
    tmp <- order(tmp, genotypeOrder=genotypeOrder(object))
    retval$genotype.freq <- retval$genotype.freq[tmp, ]

    ### from code submitted by David Duffy <davidD@qimr.edu.au>
    #
    n.typed<-sum(gf)
    correction<-2*n.typed/max(1,2*n.typed-1)
    ehet<-(1-sum(paf*paf))
    matings<- (paf %*% t(paf))^2
    uninf.mating.freq <- sum(matings)-sum(diag(matings))
    pic<- ehet - uninf.mating.freq

    retval$Hu <- correction * ehet
    retval$pic <- pic
    retval$n.typed <- n.typed
    retval$n.total <- length(object)
    retval$nallele <- nallele(object)
    #
    ###

    ## Add info on NA values
    if(any(is.na(object))){
        retval$allele.freq    <- rbind(retval$allele.freq,
                                       "NA"=c(sum(is.na(allele(object))),NA))
        retval$genotype.freq  <- rbind(retval$genotype.freq,
                                       "NA"=c(sum(is.na(object)),NA))
    }

    return(retval)
  }

print.summary.genotype  <-  function(x,...,round=2)
  { 
    if(!is.null(x$locus))
      {
        cat("\n")
        print( x$locus )
      }

    cat("\n")
    cat("Number of samples typed: ", x$n.typed,
      " (", round(100*x$n.typed/x$n.total,1), "%)\n", sep="")

    cat("\n")
    cat("Allele Frequency: (", x$nallele, " alleles)\n", sep="")
    print(round(x$allele.freq,digits=round),...)
    cat("\n")


    cat("\n")
    cat("Genotype Frequency:\n")
    print(round(x$genotype.freq,digits=round),...)
    cat("\n")

    cat("Heterozygosity (Hu)  = ",  x$Hu, "\n", sep="")
    cat("Poly. Inf. Content   = ",  x$pic, "\n", sep="")
    cat("\n")

    invisible(x)
  }
