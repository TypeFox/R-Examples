####################################
#
# Analysis of Molecular Variance
# 
# Luke Pembleton
# luke.pembleton@ecodev.vic.gov.au
#
###################################

stamppAmova <- function(dist.mat, geno, perm=100){
    
  if(class(geno)=="genlight"){  #if input file is a genlight object convert to a data.frame
    
    geno2 <- geno
    
    geno <- as.matrix(geno2) #extract genotype data from genlight object
    sample <- row.names(geno) #individual names
    pop.names <- pop(geno2) #population names
    ploidy <- ploidy(geno2) #ploidy level
    geno=geno*(1/ploidy) #convert genotype data (number of allele 2) to precentage allele frequency
    geno[is.na(geno)]=NaN
    format <- vector(length=length(geno[,1])) 
    format[1:length(geno[,1])]="genlight"
    
    
    pops <- unique(pop.names) #population names
    
    pop.num <- vector(length=length(geno[,1])) #create vector of population ID numbers
    
    for (i in 1:length(geno[,1])){
      pop.num[i]=which(pop.names[i]==pops) #assign population ID numbers to individuals
    }
    
    genoLHS <- as.data.frame(cbind(sample, pop.names, pop.num, ploidy, format))
    
    geno <- cbind(genoLHS, geno) #combine genotype data with labels to form stampp geno file
    
    geno[,2]=as.character(pop.names)
    geno[,4]=geno2@ploidy
    
    row.names(geno)=NULL
    
  }
  
  pop.names <- geno[,2]
  
  pop.names <- factor(pop.names) #updated line for compatibility with pegas 0.6
  
  temp <- environment(environment) #create temp environment
  
  assign("dist", dist.mat, envir=temp)
  assign("pop.names", pop.names, envir=temp)
  assign("perm", perm, envir=temp)
    
  res <- with(temp, amova(dist ~ pop.names, nperm=perm))
  
  rm(pop.names, perm, temp) #remove temp env
  
  return(res)
  
  
}