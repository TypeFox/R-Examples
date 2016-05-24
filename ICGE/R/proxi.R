proxi<- function(d, dx0, pert="onegroup"){
#################### Calculates the proximities #############################
# Input:
# d: distance matrix between individuals (nxn) 
# dx0: vector of length n with distancies from the specific individual to the
#      individuals of different groups.
# pert: integer vector indicating the group each individual belongs to.
#
# Output: 
# phi: vector of proximities from the specific individual to each 
#      cluster
##############################################################################


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

# We need teh squared distances
d <- d*d
dx0 <- dx0*dx0




phi <- matrix(0,k,1)

frec <- tabulate(pert) # vector of frecuencies of individuals in each population

for (pob in 1:k){
    phi[pob] <- sum(dx0[pert==pob])
    phi[pob] <- phi[pob]/frec[pob]-var[pob]
}


dimnames(phi) <- list(group=1:k, "proximity")

return(phi)

}
