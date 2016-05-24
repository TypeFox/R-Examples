INCAindex <- function(d,pert_clus){
############### Number of atypical individuals  #############################
# Input:
# d: distance matrix (n x n)
# pert_clus: partition of individuals  (n x 1)
# Output:
# Atypicals: Number of atipicals or well classified for each cluster given a 
#            partition pert_clus
# Total : Mean percentage of atypicals(well classified) for the partition 
#         pert_clus
# Ni_cluster: number of individuals in each cluster for the given partition
#############################################################################

d <- as.matrix(d)
n <- dim(d)[1]
pert_clus <- as.integer(pert_clus)
nclus<- max(pert_clus)
# populations must be named with numbers from 1 to k
if (length(tabulate(as.factor(pert_clus))) != nclus)
  stop("Partitions must be named by factors or with numbers from 1 to k.")
# 0 can not be a partitions name
if (any(pert_clus==0))
  stop("pert contains 0 named individuals.Partitions must be named by factors or with numbers from 1 to k.")

k <- nclus-1 # When there are 5 clusters, INCA is applied for 4 clusters
atypicals <- matrix(0,nclus,1)

percent_aty <- matrix(0, nclus,1)


f <- tabulate(pert_clus) #frecuencies of each cluster

   
# Calculate W for ind. in each cluster and afterwards to compare with the rest
    
 for (tt in 1:nclus){    # se va a calcular INCA para el cluster tt
     # Calculate frecuencies in each cluster to apply INCA
     ff <- matrix(0,k,1)
     aux <-0
     for (i in 1:nclus){
        if (i !=tt){
           aux <- aux+1
           ff[aux]<-f[i]
        }
     }

     
     nn <- sum(ff)  # total  individual in INCA-typicality
     nv <- f[tt] # number of individual to be testes by INCA-typicality

     pert <- matrix(0, nn,1) # partition in k clusters
     aux <- 0
     for (i in 1:n){
        if (pert_clus[i] != tt){
           aux <- aux +1
           if (pert_clus[i]<=tt){
               pert[aux] <- pert_clus[i]
           }
          if (pert_clus[i]>tt){
               pert[aux] <- pert_clus[i]-1
           }
        }
     }

      selec_xx <- pert_clus != tt
      selec_v <- pert_clus == tt
  
   
      xx_dist <- d[selec_xx,selec_xx ] 
      v_dist <- d[selec_v,selec_xx ]
  
# Verified.

       
        
     # Calculate INCA to each ind. in v_dist with respect data in xx_dist
        if (k==1){
           atypicals[tt] <- maxW_k1(xx_dist, v_dist)
           percent_aty[tt] <- atypicals[tt]/f[tt]
        }
         


       if (k>1){
           atypicals[tt] <- maxW_k(xx_dist, v_dist, pert)
           percent_aty[tt] <- atypicals[tt]/f[tt]
        }
       
        
   } # for (tt in 1:nclus)

total <- sum(percent_aty)/nclus


out <- list(well_class=atypicals, Ni_cluster=f, Total= total)
class(out) <- "incaix"

return(out)

} #end of function
