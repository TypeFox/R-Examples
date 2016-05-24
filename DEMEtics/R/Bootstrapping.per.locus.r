Bootstrapping.per.locus <- function(table2,HWEs){
  # Function used within the functions Bootstrapping.p.r and
  # Bootstrapping.CI.r to resample alleles/genotypes
                                                      if (HWEs[which(names(HWEs)==as.character((table2$locus)[1]))]==TRUE){
                          
                                                                # The confidence limits of the measure of genetic distance for the
                                                                # several loci and over all loci are determined by a thousandfold
                                                                # resampling of the alleles (for each locus
                                                                # and all populations) if the populations are in Hardy Weinberg equilibrium.
                                                                # Literature: Goudet J, Raymond M, deMeeus T, Rousset F. (1996).
                                                                # Testing differentiation in diploid populations. Genetics 144,1933-1940.
                                                                                                
                                                                            allelepool<-as.numeric(as.vector(table2$fragment.length))

                                                                                      # The alleleles that have been found at a locus in all the populations
                                                                                      # are collected in a common vector named 'allelepool'.
                                                                            
                                                                            resampled.allelepool<-sample(allelepool,length(allelepool),replace=TRUE)

                                                                                      # The alleles from one locus that were found in all populations that
                                                                                      # were sampled, are resampled with replacement.

                                                                            table2$fragment.length<-resampled.allelepool

                                                                                      # The alleles for the actual locus are replaced with alleles from
                                                                                      # the resampled.allelepool.

                                                                            table2#tab3<-rbind(tab3,tab2[[l]])

                                                                                      # The tables with the per locus resampled alleles are bound together
                                                                                      # to a single data frame.
                                                                            }else{
                                          
                                                                                     # The confidence limits of the measure of genetic distance for the
                                                                                     # several loci and over all loci are determined by a thousandfold
                                                                                     # resampling of the genotypes (for each locus and all populations) 
                                                                                     # if the populations are not in Hardy Weinberg equilibrium.
                                                                                     # Literature: Goudet J, Raymond M, deMeeus T, Rousset F. (1996).
                                                                                     # Testing differentiation in diploid populations. Genetics 144,1933-1940.
                                                                                                                                
                                                                                   table2.genotype <- split(table2$fragment.length,as.vector(table2$individual))

                                                                                              # The genotypes found for the actual locus are filtered out of
                                                                                              # table2.  They are now represented according to the frequency
                                                                                              # with which they occured in the empirical data.

                                                                                   number.genotypes <- length(table2.genotype)

                                                                                              # The number of genotypes that have to be resampled.

                                                                                   genotypepool<-as.data.frame(as.matrix(table2.genotype))[1:number.genotypes,]

                                                                                              # The genotypes that have been found for a locus in all the populations
                                                                                              # are collected in a common vector named 'allelepool'.

                                                                                   resampled.genotypepool<-sample(genotypepool,number.genotypes,replace=TRUE)
                                                                                   resampled.genotypepool <- as.numeric(as.vector(unlist(resampled.genotypepool)))

                                                                                   table2$fragment.length <- resampled.genotypepool

                                                                                              # The genotypes for the actual locus are replaced with the genotypes from
                                                                                              # the resampled genotypepool.

                                                                                 table2#  tab3<-rbind(tab3,tab2[[l]])

                                                                                              # The tables with the per locus resampled alleles are bound together
                                                                                              # to a single data frame.
          
                                                                                  }
                                                                                                                       
                                                    }
                                                      
