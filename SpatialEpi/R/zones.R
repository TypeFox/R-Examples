zones <-
function(geo, population, pop.upper.bound){

# number of areas
n <- nrow(geo)
# total population
total.pop <- sum(population)
# Interpoint distance matrix
dist <- as.matrix(dist(as.matrix(geo), upper=TRUE, diag=TRUE))


#-------------------------------------------------------------------------------
# For each area, list of closest neighbors up until pop.upper.bound of
# population is met.  We count the number of candidate zones in n.zones 
# Note:  for each county, the list of nearest neighbors form the total number of
# candidate zones
#-------------------------------------------------------------------------------
nearest.neighbors <- vector(mode="list", length=n)

n.zones <- 0
vector.cutoffs <- rep(0, n+1)

for(i in 1:n) {	
  # Sort the areas by distance, then include them one-by-one until cluster bound
  # is met
  neighbors <- order(dist[,i])
  
  # include only up until pop.upper.bound is hit
  neighbors <- neighbors[ which( 
      cumsum(population[neighbors])/total.pop <= pop.upper.bound
    )] 
  
  nearest.neighbors[[i]] <- neighbors
  
  # Update total number of zones
  n.zones <- n.zones + length(neighbors)
  
  # internal
  vector.cutoffs[i+1] <- n.zones	
}

cluster.coords <- cbind(rep(1:n, times=diff(vector.cutoffs)), 
                        unlist(nearest.neighbors))


#-------------------------------------------------------------------------------
# Output results
#-------------------------------------------------------------------------------
results <- list(
  nearest.neighbors = nearest.neighbors, 
  cluster.coords = cluster.coords, 
  dist=dist)

return(results)
}
