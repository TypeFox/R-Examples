
allele<-function(x, which=c(1,2) )
  {
    alleles.x  <- attr(x,"allele.map")
    retval  <- alleles.x[as.integer(x),which]
    attr(retval,"locus")  <- attr(x,"locus")
    attr(retval,"which")  <- which
    attr(retval,"allele.names")  <- allele.names(x)    
    #class(retval)  <- c("allele.genotype", class(retval))
    return( retval)
  }


