# Function used within allelefreq.r, in order to get the number of
# alleles of a specific locus in a specific population

calculate.anzahlk <- function(y3.1pop,unique.alleles){


number.of.alleles <- sapply(as.numeric(as.vector(unique.alleles)),function(x) length(which(as.numeric(as.vector(y3.1pop$fragment.length))==x)))

anzahlk <- cbind(unique.alleles,number.of.alleles,as.character(y3.1pop$population)[1],as.character(y3.1pop$locus)[1])
# The object 'anzahlk' is filled with the number of alleles of value
                                                                                                  # or length k, that occured in the actual population.

}
