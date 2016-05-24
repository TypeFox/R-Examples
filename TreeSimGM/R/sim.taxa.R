sim.taxa <-
function (numbsim, n, m=n,  distributionspname, distributionspparameters, distributionextname="rexp", distributionextparameters=0, symmetric = TRUE, complete=TRUE, labellivingsp="sp.", labelextinctsp="ext.", sampling=2, gsa=FALSE) { 
# numbsim is the number of simulated trees
# n is the Number of tips in sampled trees (Number of extant sampled leaves)
# m is the number of standing taxa that will exist on the first generated trees, to then be sampled for n number of tips. Case gsa=TRUE, m is equal to n.
# distributionspname is the name of the desired probability function that will be used for the speciation process (e.g. distributionspname <- "rexp"). Note that the name should contain an `r` before it, since it refers to the random number of the desired function (e.g. "rweibull", "runif")
# distributionspparameters are the parameters for the specific function desired for speciation. 
# IMPORTANT: this vector of fuction parameters must *start by the second one, since the first parameter will always be one for all the function and is added already by this function*. HINT: see the help of the desired function for more details (e.g. ?rexp) Example of parameter for a exponential distribution with lambda of one (distributionspparameters <- c(1)). Entry in the distributionparameters can be "#", # or c(#,#) in case of more characters
# distributionextname is the same as the distributionspname but for the probability of extinction (e.g. distributionextname <- "rexp") 
# distributionextparameters is the same as the distributionspparameters but for the extinction probability function. By default extinction is set to ZERO, i.e. no extinction (e.g. distributionextparameters <- c(0)). Entry in can be "#", # or c(#,#) in case of more characters
# symmetric tells which macro-evolutionary model should be used. If symmetric=TRUE the symmetric model will be used, else if FALSE, asymmetric model will be used. By default symmetric=TRUE
# complete: If complete = TRUE, the tree with the extinct and non-sampled lineages is returned. If complete = FALSE, the extinct and non-sampled lineages are suppressed. Complete=FALSE by default
# labellivingsp is the label that will be drawn on each tip surving until the present. An automatic sequential number will be added to the chosen name. By default labellivingsp="sp."
# labelextinctsp is the label that will be drawn on each extinct tip. By default labelextinctsp <- "ext."
# sampling: stochastic sampling, default
# gsa TRUE indicates that the sim.gsa.taxa will be used, the n parameter indicates the final number of species. Note that m needs to be allways bigger then n. If gsa = FALSE, there is no need of specifying n, once the final trees will be of size m
# entry in the distributionparameters can be "#", # or c(#,#) in case of more characters
	check<-gsa
	if (gsa==F && complete==T){check<-T}
	#print(paste("gsa=", gsa, " complete=", complete, " check=", check))
	mytreegsazed <- list()
while (length(mytreegsazed) < numbsim)
{
	mytree <- list()
	step <- 1
	{
	if (symmetric == TRUE) 	
	{
		for (step in 1: (numbsim) ){
			mytreenext <- mytree.symmetric.taxa(m=m, distributionspname=distributionspname, distributionspparameters=distributionspparameters, distributionextname=distributionextname, 	distributionextparameters=distributionextparameters, complete=check, labellivingsp=labellivingsp, labelextinctsp=labelextinctsp)
			mytree<- c(mytree, list(mytreenext))
			#print(paste("using symmetric [", step, "] - generated livingsp:", length(mytreenext[2]$tip.label) ))
		}	
	}
	else
	{
			for (step in 1: (numbsim) ){
			mytreenext <- mytree.asymmetric.taxa(m=m, distributionspname=distributionspname, distributionspparameters=distributionspparameters, distributionextname=distributionextname, 	distributionextparameters=distributionextparameters, complete=check, labellivingsp=labellivingsp, labelextinctsp=labelextinctsp)
			mytree<- c(mytree, list(mytreenext))
			#print(paste("using asymmetric [", step, "] - generated livingsp:", length(mytreenext[2]$tip.label) ))
	
		}
	}
	}
	#print(paste("starting gsa.taxa on", length(mytree), "trees."))
	{
	if (gsa==T)
	{
		#print("YES - gsa")
		mytreegsa <- sim.gsa.taxa(mytree, n=n, sampling=sampling, complete=complete)
	} 
	else 
	{
		#print("NO - gsa")
		mytreegsa <- mytree 
	}
	}
	
	mytreegsazed <- c(mytreegsazed, mytreegsa)
	#print(paste("finished gsa.taxa with size of ",length(mytreegsazed), "trees"))
	}
mytreegsazeds <- sample(mytreegsazed, numbsim)
mytreegsazed <- mytreegsazeds
return(mytreegsazed)	
}
