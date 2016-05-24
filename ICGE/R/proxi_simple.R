proxi_simple<- function(dx0, var, pert){
################### Proximity function (faster version)########################
# Input:  
# dx0: vector of dist. from one individual to the ind. of different clusters.
# var: geometric variabilities.
# pert: integer vector indicating the group each individual belongs to.
# Output: 
# phi: vector of proximities from one individual to each cluster.

dx0 <- dx0*dx0

n <- length(pert)
k<- max(pert)


phi <- matrix(0,k,1)

frec <- tabulate(pert) # vector of frecuencies of individuals in each population


for (pob in 1:k){
    phi[pob] <- sum(dx0[pert==pob])
    phi[pob] <- phi[pob]/frec[pob]-var[pob]
}


return(phi)

}
