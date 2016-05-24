# The R code for including the graphlet count changes into the statnet package
InitErgmTerm.grorbitFactor <- function(nw, arglist, ...){
	a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite = NULL,
						varnames = c("covar" , "grorbits", "base"), 
						vartypes = c("character" , "numeric", "numeric"), 
						defaultvalues = list(NULL, c(0:72) , 1),
						required = c(TRUE, FALSE, FALSE))
	
	
	# Check the number of attributes provided for covar
	if (length(a$covar) > 1){
		print("Only one covariate can be processed at a time (at least for the moment)")
		return(NULL)
	}
	
	# Check the correctness of the provided graphlet orbit IDs 
	gorList <- c()

	for (x in a$grorbits){
		if (x > 72 || x < 0){
			print(paste("WARNING: Wrong graphlet orbit ID provided: " , x))
			print("The wrong graphlet orbit ID is ignored !")
		}
		else {
			gorList <- c(gorList , x)
		}
	}
	
	if (length(gorList) == 0){
		print("No valid graphlet orbit ID's are provided for grorbitFactor!")
		return(NULL)
	}
	
	# Identify the parameters for the attribute term
	attributes <- get.node.attr(nw, a$covar)
	
	u <- sort(unique(attributes))
	if (!is.null(a$base) && !identical(a$base,0)) {
		u <- u[-a$base]
		if (length(u)==0) {
			print("Warning:  grorbitFactor term deleted because it contributes no statistics")
			return(NULL)
		}
	}
	
	# Recode the attribute values to numeric
	nodecov <- match(attributes,u,nomatch=length(u)+1)
	
	# Define the coefficient names when the 0 statistics terms are dropped
	coef.names = c()	
	for (i in 1:length(gorList)){
		coef.names <- append(coef.names, paste("grorbitFactor.orb_" , gorList[i], ".attr_" , u , sep = ""))
	}
		
	# Identify the coefficient names 	
	list(name="grorbitFactor",
		 coef.names=coef.names,
		 inputs=c(length(gorList), gorList, length(u), network.size(nw), nodecov),
		 pkgname="ergm.graphlets",
		 dependence=TRUE, 
		 minval = 0)
}