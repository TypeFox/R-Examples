# Function called within the funciton Ht.value.calculations.r
Mean.allelefrequency <- function(Table3.1allele){
Mean <- mean(as.numeric(as.vector(Table3.1allele$proportion)))
# The mean allelefrequency over all populations is calculated for the
# actual allele.
invisible(Mean)
}
