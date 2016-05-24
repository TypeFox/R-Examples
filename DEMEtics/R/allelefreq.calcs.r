# Function called within allelefreq.r in order to calculate
# allelefrequencies for each locus and population.
allelefreq.calcs <- function(y2){
y3 <- split(y2,y2$population)

unique.alleles <- unique(y2$fragment.length) # unique alleles for the actual locus
unique.alleles <- unique.alleles[unique.alleles!=0] # 0 is not an allele but indicates non-amplification, therefore it is removed from the vector containing allele lengths.

anzahlk1 <- lapply(y3, function(x) calculate.anzahlk(x,unique.alleles))
anzahlk <- do.call(rbind,anzahlk1)
anzahlk=as.data.frame(anzahlk)
colnames(anzahlk) <- c("allele","number","population","locus")
# The object 'anzahlk' is set to a data frame format and the columns
                                      # of this object are named.
anzahlk.pop <- split(anzahlk,anzahlk$population)

percents <- lapply(anzahlk.pop,pop.allele.proportions)

proportion=as.numeric(as.vector(do.call(c,percents)))

allelefrequency <- cbind(anzahlk,proportion)
allelefrequency                                                                                                  
}
