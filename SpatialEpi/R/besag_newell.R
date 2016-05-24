besag_newell <-
function(geo, population, cases, expected.cases=NULL, k, alpha.level){

#-------------------------------------------------------------------------------
# Initialization 
#-------------------------------------------------------------------------------
# If no expected.cases provided, set them if there are no expected counts
if(is.null(expected.cases)){
	p <- sum(cases)/sum(population)
	expected.cases <- population*p
}

# geographical information computation
geo.results <- zones(geo, population, 1)
nearest.neighbors <- geo.results$nearest.neighbors
distance <- geo.results$dist
n.zones <- length(unlist(nearest.neighbors))


#-------------------------------------------------------------------------------
# Observed statistic computation
#-------------------------------------------------------------------------------
results <- besag_newell_internal(cases, expected.cases, nearest.neighbors, 
                                 n.zones, k)

# observed p.values for each areas
p.values <- results$observed.p.values	
# observed number of neighbors needed to observe k cases
m.values <- results$observed.m.values	
# actual observed number of cases
k.values <- results$observed.k.values	

# pick out areas that were significant and order them by p-value
signif.indices <- order(p.values)[1:sum(p.values <= alpha.level)]

# order remaining values
signif.p.values <- p.values[signif.indices]
signif.m.values <- m.values[signif.indices]
signif.k.values <- k.values[signif.indices]


# Create object to output
# If none are significant, return NULL
if(length(signif.indices) == 0){
	clusters <- NULL
} else {
	clusters <- vector("list", length=length(signif.indices))

	for( i in 1:length(clusters) ){	
		# find areas included in cluster
		cluster <- order(distance[signif.indices[i],])[1:signif.m.values[i]]
		
		new.cluster <- list(
			location.IDs.included = cluster,
			population = sum(population[cluster]),
			number.of.cases = sum(cases[cluster]),
			expected.cases = sum(expected.cases[cluster]),
			SMR = sum(cases[cluster])/sum(expected.cases[cluster]),
			p.value = signif.p.values[i]	
		)
		clusters[[i]] <- new.cluster
	}
}


#-------------------------------------------------------------------------------
# Output results
#-------------------------------------------------------------------------------
results <- list(
	clusters=clusters,
	p.values=p.values,
	m.values=m.values,
	observed.k.values=k.values
)	
return(results)
}
