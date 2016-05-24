# The R code for including the graphlet count changes into the statnet package
InitErgmTerm.graphletCount <- function(nw, arglist, ...){
	a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite = NULL,
						varnames = c("graphlets"), 
						vartypes = c("numeric"), 
						defaultvalues = list(0:29),
						required = c(FALSE))
	
	grList <- c()
	
	for (x in a$graphlets){
		if (x > 29 || x < 0){
			print(paste("WARNING: Wrong graphlet ID provided: " , x))
			print("The wrong parameter is ignored")
		}
		else {
			grList <- c(grList , x)
		}
	}
	
	if(length(grList) == 0){
		print("No valid graphlet ID is provided! Evaluating all graphlets...")
        grList <- 0:29
	}
	
	coef.names <- paste("graphlet." , grList, ".Count", sep = "")
	
	list(name="graphlets",
		 coef.names=coef.names,
		 inputs=c(length(grList),grList),
		 pkgname="ergm.graphlets",
		 dependence=TRUE, 
		 minval = 0)
}