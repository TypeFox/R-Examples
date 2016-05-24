kulldorff <-
function(geo, cases, population, expected.cases=NULL, pop.upper.bound, 
         n.simulations, alpha.level, plot=TRUE){

#-------------------------------------------------------------------------------
# Initialization 
#-------------------------------------------------------------------------------
# Determine likelihood type: binomial or poisson
if(is.null(expected.cases)){
	type <- "binomial"
	denominator <- population
	expected.cases <- sum(cases) * (denominator/sum(denominator))
# poisson case	
}else{
	type <- "poisson"
	denominator <- expected.cases
}

# Get geographic information
geo.results <- zones(geo, population, pop.upper.bound)
nearest.neighbors <- geo.results$nearest.neighbors
cluster.coords <- geo.results$cluster.coords
n.zones <- nrow(cluster.coords)


#-------------------------------------------------------------------------------
# Observed statistic computation
#-------------------------------------------------------------------------------
lkhd <- computeAllLogLkhd(cases, denominator, nearest.neighbors, n.zones, type)

# Get areas included in most likely cluster
cluster.index <- which.max(lkhd)

# cluster center and radial area
center <- cluster.coords[cluster.index,1]
end <- cluster.coords[cluster.index,2]

# list of all areas included in cluster	
cluster <- nearest.neighbors[[center]]
cluster <- cluster[1:which(cluster == end)]


#-------------------------------------------------------------------------------
# Compute Monte Carlo randomized p-value
#-------------------------------------------------------------------------------
# Simulate cases under null hypothesis of no area effects i.e. conditioned on E
perm <- rmultinom(n.simulations, round(sum(cases)), prob=denominator)

# Compute simulated lambda's:  max log-lkhd in region
sim.lambda <- kulldorffMC(perm, denominator, nearest.neighbors, n.zones, type)

# Compute Monte Carlo p-value
combined.lambda <- c(sim.lambda, max(lkhd))
p.value <- 1-mean(combined.lambda < max(lkhd))

# Plot histogram
if(plot){
	hist(combined.lambda, 
       main="Monte Carlo Distribution of Lambda",
       xlab=expression(log(lambda)))
	abline(v=max(lkhd), col="red")
	legend("top",
		c(paste("Obs. log(Lambda) = ",round(max(lkhd),3),sep=""), 
			paste("p-value = ", round(p.value,log10(n.simulations + 1)),sep="")),
		lty=c(1, 1), 
		col=c("red","white"), 
		bty="n"
	)
}


#-------------------------------------------------------------------------------
# Create Most Likely Cluster Object
#-------------------------------------------------------------------------------
most.likely.cluster = list(
	location.IDs.included = cluster,
	population = sum(population[cluster]),	
	number.of.cases = sum(cases[cluster]),
	expected.cases = sum(expected.cases[cluster]),
	SMR = sum(cases[cluster])/sum(expected.cases[cluster]),
	log.likelihood.ratio = lkhd[cluster.index],
	monte.carlo.rank = sum(combined.lambda >= lkhd[cluster.index]),
	p.value = p.value
)


#-------------------------------------------------------------------------------
# Investigate Secondary Clusters
#-------------------------------------------------------------------------------
# running list of areas already covered to test overlap
current.cluster <- cluster
secondary.clusters <- NULL

# go through log-likelihoods in decreasing order, skipping largest which 
# corresponds to most likely cluster
indices <- order(lkhd, decreasing=TRUE)

for(i in 2:length(indices)){
	#---------------------------------------------------
	# Get areas included in new potential cluster
	new.cluster.index <- indices[i]
	new.center <- cluster.coords[new.cluster.index, 1]
	new.end <- cluster.coords[new.cluster.index, 2]
	
	new.cluster <- nearest.neighbors[[new.center]]
	new.cluster <- new.cluster[1:which(new.cluster == new.end)]	
	#---------------------------------------------------
	# if there is no overlap between existing clusters and the new potential cluster
	if(length(intersect(new.cluster,current.cluster)) == 0){		
		new.secondary.cluster <- list(
			location.IDs.included = new.cluster,
			population = sum(population[new.cluster]),
			number.of.cases = sum(cases[new.cluster]),
			expected.cases = sum(expected.cases[new.cluster]),
			SMR = sum(cases[new.cluster])/sum(expected.cases[new.cluster]),
			log.likelihood.ratio = lkhd[new.cluster.index],
			monte.carlo.rank = sum(combined.lambda >= lkhd[new.cluster.index]),
			p.value =  1 - mean(combined.lambda < lkhd[new.cluster.index])
		)

		# If no longer significant, end loop.
		if(new.secondary.cluster$p.value > alpha.level){
			break
		}
		# Otherwise add to existing clusters
		secondary.clusters <- c(secondary.clusters,list(new.secondary.cluster))
		current.cluster <- unique(c(current.cluster, new.cluster))
	}		
}


#-------------------------------------------------------------------------------
# Output results
#-------------------------------------------------------------------------------
results <- list(
	most.likely.cluster=most.likely.cluster,
	secondary.clusters=secondary.clusters,
	type = type,
	log.lkhd=lkhd,
	simulated.log.lkhd=sim.lambda
)
return(results)
}
