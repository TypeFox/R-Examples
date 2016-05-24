# Function called within H.out.r

Hj.values.calculation <- function(allelefrequency2.onelocus){
allelefrequency3<-split(allelefrequency2.onelocus,allelefrequency2.onelocus$population)
# For the actual locus, the table allelefrequency2 is splitted according
# to the several populations.
  
Hj.values.onelocus1 <- lapply(allelefrequency3,Hj.values.one.population)
Hj.values.onelocus <- do.call(rbind,Hj.values.onelocus1)
invisible(Hj.values.onelocus)

}
