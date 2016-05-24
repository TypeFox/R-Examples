# Function used in Hardy.Weinberg.r to calculate the estimated expected and empirically obtained number of certain genotypes

Frequency.genotypes.calculation <- function(seq.pair.comb,tab2.pop.genotypes,tab2.pop.alleles,allele.one,allele.two,tab2.pop.number.individuals){

first.allele <- as.numeric(as.vector(tab2.pop.alleles[allele.one[seq.pair.comb],1]))
# The first allele of the actual genotype that is searched for
second.allele <- as.numeric(as.vector(tab2.pop.alleles[allele.two[seq.pair.comb],1]))
# The second allele of the actual genotype that is searched for.


estimated.frequency.genotypes <- ifelse(first.allele==second.allele,(tab2.pop.alleles[allele.one[seq.pair.comb],2])^2,2*(tab2.pop.alleles[allele.one[seq.pair.comb],2])*tab2.pop.alleles[allele.two[seq.pair.comb],2])

tab2.pop.seq <- seq(1:tab2.pop.number.individuals)

this.genotype1 <- lapply(tab2.pop.seq, function(x) this.genotype.calc(x,tab2.pop.genotypes,first.allele,second.allele))
this.genotype <- do.call(c,this.genotype1)

                                                                                             
empirical.number.genotypes <- length(which(this.genotype==TRUE))
                                                                  
# The number of occurences of the several genotypes are combined in
# this vector.

out <- cbind(estimated.frequency.genotypes,empirical.number.genotypes)
invisible(out)
}
