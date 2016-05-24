# This function is used within the function
# Frequency.genotypes.calculation.r.  It tests for the presence of a
# specific allele combination and puts out a TRUE value if the
# combination is present, otherwise it returns FALSE.
this.genotype.calc <- function(tab2.pop.seq,tab2.pop.genotypes,first.allele,second.allele){
# This will calculate the number of individuals that possess the
# actual genotype and place the results in a vector.
  
# Now, in every individual, the actual allele combination is searched.
                              
occurrence <- as.numeric(as.vector(tab2.pop.genotypes[[tab2.pop.seq]]$fragment.length[1]))==first.allele&&as.numeric(as.vector(tab2.pop.genotypes[[tab2.pop.seq]]$fragment.length[2]))==second.allele || as.numeric(as.vector(tab2.pop.genotypes[[tab2.pop.seq]]$fragment.length[2]))==first.allele&&as.numeric(as.vector(tab2.pop.genotypes[[tab2.pop.seq]]$fragment.length[1]))==second.allele

 # If the actual allele combination is represented in the actual
 # individual, a TRUE is returned, otherwise a FALSE.
          
}
