# Function called within Bootstrapping.CI.r and in Bootstrapping.p.r
locus.means.calc <- function(repetition.val,tab2,HWEs,x){

table3 <- lapply(tab2,function(y) Bootstrapping.per.locus(y,HWEs))
                            len <- length(table3)
                            table3 <- as.data.frame(do.call(rbind, table3))                                                                                                                       
                            allelefreq(table3)
                            
                                      # The table that contains the allelefrequencies and the table that
                                      # lists the sample sizes are placed in the R workspace in the 
                                      # object List, but also separately in the object allelefrequency
                                      # and the object sample.sizes by this function.                              
                            
                            calc(DEMEtics.env$allelefrequency,DEMEtics.env$sample.sizes,x)
                            
                                      # The Dest values for the several loci and over all loci are
                                      # calculated for all the populations that have been examined.
                                      # The results are available from the object 'D.values'. 
                            
                            locus <-DEMEtics.env$values[[1]]
                            means <-DEMEtics.env$values[[2]]
out <- list(locus,means)
invisible(out)
  
}
