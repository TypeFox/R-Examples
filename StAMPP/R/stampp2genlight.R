#####################################################
#
# Convert StAMPP genotype object to Genlight object
# 
# Luke Pembleton
# luke.pembleton@ecodev.vic.gov.au
# 
######################################################

stampp2genlight <-
function (geno, pop=TRUE){
    
                    
      data <- geno[,-(1:5)]   #matrix of allele frequencies
      ind <- geno[,1] #individual ids
      ploidy.levels <- geno[,4] #ploidy
      pop.names <- geno[,2] #population ids
      
      data=data*ploidy.levels #convert percentage of allele A to number of allele A based on ploidy level
      
      data <- new("genlight", data, parallel=FALSE) #convert genotype data to genlight object
      indNames(data)=ind #add individual ids to genlight object
      ploidy(data)=ploidy.levels #add ploidy levels to genlight object
      
      if(pop==TRUE){
        
        pop(data)=pop.names #if population ids are present, add to genlight object
        
      }
        
    
    return(data)
    
  }
