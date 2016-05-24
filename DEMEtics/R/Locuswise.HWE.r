Locuswise.HWE <- function(table.of.list){
  # Function used within Bootstrapping.p.r and in Bootstrapping.CI.r
  # to identify if all populations are for the respective locus in HWE.
  
                         HWE <- Hardy.Weinberg(table.of.list,1)
                          
                                    # It is tested, if all populations are in Hardy Weinberg equilibrium
                                    # for the actual locus.
                                    # The result is either HWE=TRUE or HWE=FALSE.                 

                           #HWEs[l] <- HWE  
                           
                                    # The results for the several loci are combined in a single vector.
                                    
                           if (HWE==TRUE){
                                          cat("\n","All of these populations are in Hardy Weinberg Equilibrium with regard to the locus: ",as.character((table.of.list$locus[1])),"\n",sep="")
                                          cat("Therefore, alleles are permuted among these populations for this locus.","\n")                                      
                                          }else{
                                          cat("\n","Not all of these populations are in Hardy Weinberg Equilibrium with regard to the locus: ",as.character((table.of.list$locus)[1]),"\n",sep="")
                                          cat("Therefore, genotypes are permuted among these populations for this locus.","\n","\n")                                      
                                          }

HWE                                                    # User information about the permutation method and its reasons.                                          
                                                                     
                       }
