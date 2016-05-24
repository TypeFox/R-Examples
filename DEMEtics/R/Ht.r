Ht <- function(Table){

#Variables
#------------------------------------------------------------------------------------------------------------------------------
# Input:
#           Table <- Dest.calc(), Gst.est.calc();
#------------------------------------------------------------------------------------------------------------------------------
  

          # function to calculate the total heterozygosity.
    
          # The argument 'Table' has to be of the following format:

          #   allele number population   locus  proportion

          #       4      1     Borkum  LoPi89 0.1250000
          #       5      1     Borkum  LoPi89 0.1250000
          #       6      2     Borkum  LoPi89 0.2500000
          #       7      0     Borkum  LoPi89 0.0000000
          #       8      3     Borkum  LoPi89 0.3750000
          #       9      1     Borkum  LoPi89 0.1250000
          #       4      1   Langeoog  LoPi89 0.1666667
          #       5      1   Langeoog  LoPi89 0.1666667
          #       6      1   Langeoog  LoPi89 0.1666667
          #       7      2   Langeoog  LoPi89 0.3333333
          #       8      1   Langeoog  LoPi89 0.1666667
          #       9      0   Langeoog  LoPi89 0.0000000
          #       4      0 Wangerooge  LoPi89 0.0000000
          #       5      1 Wangerooge  LoPi89 0.1250000
          #       6      0 Wangerooge  LoPi89 0.0000000
          #       7      4 Wangerooge  LoPi89 0.5000000
          #       8      2 Wangerooge  LoPi89 0.2500000
          #       9      1 Wangerooge  LoPi89 0.1250000


Table2 <- split(Table,Table$locus)
Ht.one.locus <- lapply(Table2,Ht.value.calculations)
Ht.values1 <- do.call(rbind,Ht.one.locus)
Ht.values <- cbind(Ht.values1,names(Table2))
# The Ht values for the several loci are combined with a column that
# contains the names of the loci, the Ht-value belongs to.
      
Ht.values <- as.data.frame(Ht.values)
colnames(Ht.values) <- c("Ht.value","locus")
invisible(Ht.values)
  
}
