deltas <-function(d,pert="onegroup"){
############ Calculates distances between groups ######################
# Input:
# d: distance matrix between individuals (nxn) 
# pert: integer vector indicating the group each individual belongs to.
# Output:
# delta: kxk matrix with distances between groups
########################################################################

d <- as.matrix(d)
n<-dim(d)[1]
if (pert[1]=="onegroup"){pert <- rep(1,n)}
pert <- as.integer(pert)
k<-max(pert)
# populations must be named with numbers from 1 to k
if (length(tabulate(as.factor(pert))) != k)
  stop("Partitions must be named by factors or with numbers from 1 to k.")
# 0 can not be a partitions name
if (any(pert==0))
  stop("pert contains 0 named individuals.Partitions must be named by factors or with numbers from 1 to k.")

# We need the geometrical variabilities
var <- vgeo(d,pert)

# We need the squared distances
d <- d*d


delta <- matrix(0, k,k) #matrix of distances between groups

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



dimnames(delta) <- list(group=(1:k), 1:k)

return(delta)
}

