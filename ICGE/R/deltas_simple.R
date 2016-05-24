deltas_simple <-function(d,var,pert){
############ Distances between groups (faster version) ##########################
# Input: 		
# d: distance matrix(nxn) 
# var: vector of geometric variabilities
# pert:  integer vector indicating the group each individual belongs to.
# Output:
# delta: kxk matrix with distances between populations
###############################################################################

d <- d*d

n<-dim(d)[1]
k<-max(pert)
delta <- matrix(0, k,k)

frec <- tabulate(pert) # vector of frecuencies of individuals in each population


for (pob1 in 1:(k-1)){
  for (pob2 in (pob1+1):k){
    delta[pob1, pob2] <- sum(d[pert==pob1, pert==pob2])
  }
}


for (pob1 in 1:(k-1)){
    for (pob2 in (pob1+1):k){
        delta[pob1, pob2] <- delta[pob1, pob2]/(frec[pob1]*frec[pob2])-var[pob1]-var[pob2]
        delta[pob2, pob1] <- delta[pob1, pob2]
    }
}


return(delta)
}

