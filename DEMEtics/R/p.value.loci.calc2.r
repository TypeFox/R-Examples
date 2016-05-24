# Function called within pair.pops.r

p.value.loci.calc2 <- function(loci2.1locus,values2.1locus){

bootstrapped.values.loci=as.numeric(as.vector(loci2.1locus[,1]))

# The bootstrapping values (D or Gst values) that have been obtained
# for the several loci are assigned to a vector.
          
empirical.value.loci=as.numeric(as.vector(values2.1locus[,1]))
                          
# The empirical D or Gst values for the several loci are assigned to a
# vector.                                              
                                                              
if (is.nan(empirical.value.loci)){
p.value <- NA
} else 
p.value <- p.val(empirical.value.loci,bootstrapped.values.loci)
                                                              
                                                                                                  # This function gives the p-value for the actual empirical value
# in the object 'p.value'.
                                                       
# In the case that no bootstrapped values had been obtained
# when populations that are compared contain all
# the same allele at one locus (allelefrequency 100%), a Gst or D value
# can't be calculated; (0/0) is calculated.
# The empirical value is also NA in this case and then the p.value
# is also given as NA. Otherwise a p.value is calculated from the
# bootstrapped values.
p.values.loci <- cbind(round(p.value,4),as.character((loci2.1locus$locus)[1]))                                                       
invisible(p.values.loci)
                                                            
}
