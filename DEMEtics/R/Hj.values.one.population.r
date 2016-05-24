# Function called within Hj.values.calculation.r

                                          
Hj.values.one.population <- function(allelefrequency3.onepop){

Hj.one.population<-Hj(as.numeric(as.vector(allelefrequency3.onepop$proportion)))
# The Hj-value for the actual locus and population is calculated.
                                                                  
Hj.values.onelocus<-cbind(as.character((allelefrequency3.onepop$locus)[1]),as.character((allelefrequency3.onepop$population)[1]),Hj.one.population)
                              
# The Hj-values are combined with a column of the names of the actual population
# and a column of the names of the actual locus.

invisible(Hj.values.onelocus)

}
