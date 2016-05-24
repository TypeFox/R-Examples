# The R code for including the graphlet count changes into the statnet package
InitErgmTerm.grorbitDist <- function(nw, arglist, ...){
	a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite = NULL,
						varnames = c("grorbits" , "values"), 
						vartypes = c("numeric" , "numeric"), 
						defaultvalues = list(c(0:14), NULL),
						required = c(TRUE, TRUE))
	
	# Check the correctness of the provided orbits
	gorList <- c()
	
	for (x in a$grorbits){
		if (x > 14 || x < 0){
			print(paste("WARNING: Wrong graphlet orbit ID provided to grorbitDist: " , x))
			print("The wrong graphlet orbit ID is ignored !")
		}
		else {
			gorList <- c(gorList , x)
		}
	}
	
	if (length(gorList) == 0){
		print("You should provide some orbit id's to use grorbitDist term")
		return(NULL)
	}
	
	evalVals <- a$values
	
	# Count everything in the graph once to get the distribution
	consideredOrbit <- 15	
	counts <- .C("ncount", as.integer(as.matrix(nw)) , as.integer(network.size(nw)) , as.integer(network.edgecount(nw)) , out=integer(network.size(nw) * consideredOrbit))
            
	# These variables keep track of the MCMC process
	originalCounts <- counts$out
	
	# Identify the coefficient names
	coef.names <- c()	
	
	for (i in 1:length(gorList)){
		coef.names <- append(coef.names , paste("grorbitDist.orbit_" , gorList[i] , ".val_", evalVals , sep = ""))
	}
    
	
	# Initialize the empty network statistics
	emptynwstats <- integer(length(gorList) * length(evalVals))
	
	for (x in 0:(length(gorList)-1)){
		for (y in 1:(length(evalVals))){
			if (evalVals[y] == 0) {
				emptynwstats[x * length(evalVals) + y] <- network.size(nw)
			}
		}
	}
    
	inputs <- c(length(gorList),gorList, length(evalVals), a$values, 
                    network.edgecount(nw), originalCounts)

	list(name="grorbitDist",
		 coef.names=coef.names,
		 inputs=inputs,
		 pkgname="ergm.graphlets",
		 emptynwstats=emptynwstats,
		 dependence=TRUE, 
		 minval = 0)	
}
