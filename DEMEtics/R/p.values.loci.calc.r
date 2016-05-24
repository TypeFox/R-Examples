# Funciton called within over.all.pops.r

p.values.loci.calc <- function(loci2.1locus,values2.1locus){
          
bootstrapped.values.loci=as.numeric(as.vector(loci2.1locus[,1]))

# The bootstrapping values (D or Gst values) that have been obtained
# for the several loci are assigned to a vector.
          
empirical.value.loci=as.numeric(as.vector(values2.1locus[,1]))
                          
# The empirical D or Gst values for the several loci are assigned to a
# vector.

p.value <- p.val(empirical.value.loci,bootstrapped.values.loci)

#This function returns the p-value for the actual empirical D or Gst value
# in the object 'p.value'.

          
p.values.loci=round(p.value,4)
                          
# The p-values (rounded up to 4 decimal places) can't be more exact, when the bootstrapping is based on thousand resamplings.    
}
