#Function called within pair.pops.r

v.loci.calc <- function(allelefrequency.pair2.1locus,values2.1locus){
allelefrequency.pair3=split(allelefrequency.pair2.1locus,as.character(allelefrequency.pair2.1locus$population))
                                                    
# The table 'allelefrequency.pair2 is split according to the
# populations that are actually compared It is defined again for the
# case that for one locus, data for one of the populations would be
# lacking.
                                                              
v.loci=cbind(as.numeric(as.vector(values2.1locus[,1])),as.character((allelefrequency.pair2.1locus$locus)[1]),names(allelefrequency.pair3)[1],names(allelefrequency.pair3)[2])
# In this vector,the Gst or D value for
# the actual locus, the actual locus, population one and population
# two that are actually compared with one another are combined.
invisible(v.loci)                                                                                                            
}
   
