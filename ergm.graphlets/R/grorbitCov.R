# The R code for including the graphlet count changes into the statnet package
InitErgmTerm.grorbitCov <- function(nw, arglist, ...){
	a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite = NULL,
						varnames = c("covar" , "grorbits"), 
						vartypes = c("character" , "numeric"), 
						defaultvalues = list(NULL, c(0:72)),
						required = c(TRUE, FALSE))
	
	if (length(a$covar) > 1){
		print("Only one covariate can be processed at a time")
		return(NULL)
	}
	
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
		print("No valid graphlet orbit ID's are provided !")
		return(NULL)
	}
		
	coef.names <- paste("grorbitCov" , gorList, "." , a$covar, sep = "")
	attributes <- get.node.attr(nw, a$covar, "grorbitCov", numeric=TRUE)
	
	list(name="grorbitCov",
		 coef.names=coef.names,
		 inputs=c(length(gorList), gorList, network.size(nw), attributes),
		 pkgname="ergm.graphlets",
		 dependence=TRUE, 
		 minval = 0)
}