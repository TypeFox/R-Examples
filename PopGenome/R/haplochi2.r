#Cai JJ (2008) PGEToolbox: A Matlab toolbox for population genetics and evolution
#Journal of Heredity Jul-Aug;99(4):438-40. doi:10.1093/jhered/esm127
#modified

haplochi2 <- function(populations,nh,sfreqh){ # Chi2Test

npops <- length(populations)

# Get sample size of each population
n <- vector(,npops)
for(xx in 1:npops){
   n[xx] <- length(populations[[xx]])
}

# Number of haplotypes
 K              <- dim(sfreqh)[2]
 Sumnij              <- colSums(sfreqh)
 p              <- Sumnij/K
 
 return(list(K=K,Sumnij=Sumnij,p=p))
 }
