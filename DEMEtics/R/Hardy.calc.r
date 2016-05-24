# Function used within Hardy.Weinberg.r in order to calculate the p-value for populations being in HWE for a certain locus.
Hardy.calc <- function(tab2.pop.1pop){
tab2.pop.individual <- as.character(tab2.pop.1pop$individual)
# The individuals occuring in the actual population.
tab2.pop.alleles <- as.data.frame.table(table(as.numeric(as.vector(tab2.pop.1pop$fragment.length))))
# The number of the several alleles of the actual locus occuring in the 
                                    # actual population. 
number.alleles <- tab2.pop.alleles[,2]      
all.alleles <- sum(tab2.pop.alleles[,2])
frequency.alleles <- number.alleles/all.alleles
                          
# The frequency of the several alleles is calculated.

tab2.pop.alleles[,2] <- frequency.alleles
                          
# The numbers of the different alleles are replaced by their 
# frequencies with which they occur in the actual population.

tab2.pop.number.individuals <- length(tab2.pop.individual)/2

# The number of individuals occuring in the actual population.

tab2.pop.genotypes <- split(as.data.frame(as.matrix(tab2.pop.1pop)),as.data.frame(as.matrix(tab2.pop.individual)))
                          
# The genotypes occuring for the actual locus and population.  For
                          # each population it has to be found out
                          # seperately, if it is in Hardy Weinberg
                          # equilibrium for the actual locus.
                          
tab2.pop.number.alleles <- length(tab2.pop.alleles[,1])
                          
# The number of alleles for the actual locus and population.
          
pairwise.combinations <- sum(seq(1:(tab2.pop.number.alleles)))
                          
# Depending on the number of alleles, the number of possible
# combinations of alleles is calculated.
          
# Construction for pairwise combination:

repetition <- seq(tab2.pop.number.alleles,1)
                          
allele.one1 <- lapply(repetition, function(x) allele.one <- rep(max(repetition)+1-x,x))
allele.one <- do.call(c,allele.one1)

# With these commands, the positions of the first allele in the
# table tab2.pop.alleles that will be combined with another
# allele, are created.
                          
allele.two1 <- lapply(repetition, function(x) allele.two <- seq(max(repetition)+1-x,max(repetition)))
allele.two <- do.call(c,allele.two1)
# The vector of the positions of the alleles in the table tab2.pop.alleles,
# with which the first chosen alleles in the vector 'allele.one'
# will be compared with, is created.

empirical.number.genotypes <- numeric(0)
                          
# This vector will be filled with the number of occurences of the
# several genotypes for the actual locus in the actual population.
                                    
estimated.frequency.genotypes <-  numeric(0)
# This vector will be filled with the estimated frequency of genotypes
# that would be expected if the population would be in HWE for the
# actual locus.
seq.pair.comb <- seq(1,pairwise.combinations)

counted.genotypes1 <-  lapply(seq.pair.comb, function(x) Frequency.genotypes.calculation(x,tab2.pop.genotypes,tab2.pop.alleles,allele.one,allele.two,tab2.pop.number.individuals))

counted.genotypes <- do.call(rbind,counted.genotypes1)
estimated.frequency.genotypes <- as.numeric(as.vector(counted.genotypes[,1]))
empirical.number.genotypes <- as.numeric(as.vector(counted.genotypes[,2]))
                                                   
estimated.number.genotypes <- tab2.pop.number.individuals*estimated.frequency.genotypes
                           
 # The estimated number of genotypes, if the population were in HWE
                                    # for the actual locus, is
                                    # calculated.
          
 if (length(estimated.number.genotypes)==1 || all(estimated.number.genotypes==estimated.number.genotypes[1]) ||  all(empirical.number.genotypes==empirical.number.genotypes[1])) {
                                                                      p.value <- 1
                                                                      } else {
                                                                      p.value <- chisq.test(estimated.number.genotypes,empirical.number.genotypes,simulate.p.value=TRUE,B=10000)$p.value
                                                                      }                                                                      
 # If a locus is made up of just one single allele, a chisquare test
 # can't be carried out. In this case, with p.value <- 1, the actual
 # population is set to be in Hardy Weinberg equilibrium for this
 # locus.
          
p.value <- ifelse (is.na(p.value), 0, p.value)       
                          
# In the case that a Chisquare test can't be carried out, the randomization
# of alleles is not justified. Therefore, in this case, the p.value is set
# to 0 (not in HW-Equilibrium).                                                                              
invisible(p.value)          
            

                            
}
