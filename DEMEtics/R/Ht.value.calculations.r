# Function called within Ht.r
Ht.value.calculations <- function(Table2.1locus){
Table3 <- split(Table2.1locus,as.numeric(as.vector(Table2.1locus$allele)))
# The allelefrequency-values are splitted according to the several
# alleles that were found in this population
Means <- as.numeric(as.vector(sapply(Table3,Mean.allelefrequency)))
Homozygotes <- Means^2
Ht.one.locus <- 1-sum(Homozygotes)
invisible(Ht.one.locus)
}
