vgeo<-function(d,pert="onegroup"){
############ Calculates the geometric variabilities ##################
# Input:
# d: distance matrix between n individuals
# pert: integer vector indicating the group each individual belongs 
#
# Output: vector of geometrical variabilities of each group
#######################################################################

d <- as.matrix(d)
n <- dim(d)[1] # there are n elements
if(pert[1]== "onegroup") {pert <- rep(1,n)}
pert <- as.integer(pert)
k <- max(pert) # there are k populations
# populations must be named with numbers from 1 to k
if (length(tabulate(as.factor(pert))) != k)
  stop("Partitions must be named by factors or with numbers from 1 to k.")
# 0 can not be a partitions name
if (any(pert==0))
  stop("pert contains 0 named individuals.Partitions must be named by factors or with numbers from 1 to k.")

# We need squared distances
d <- d*d


var<-matrix(0,k,1) #vector of geometric variabilities 


frec <- tabulate(pert)

for (pob in 1:k){
     var[pob] <- sum(d[pert == pob , pert == pob])
}


var <- var/(2*frec*frec)

  
dimnames(var) <- list(group=(1:k), "geom. variability")


return(var)  

}
