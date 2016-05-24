sim.age <-
function (age, numbsim, distributionspname, distributionspparameters, distributionextname="rexp", distributionextparameters=0, symmetric = TRUE, complete=TRUE, labellivingsp="sp.", labelextinctsp="ext.") { 
# age is the total age until each tree will be simulated (e.g. age <- 3). Time since origin / most recent common ancestor
# numbsim is the number of simulated trees
# distributionspname is the name of the desired probability function that will be used for the speciation process (e.g. distributionspname <- "rexp"). Note that the name should contain an `r` before it, since it refers to the randon number of the desired function (e.g. "rweibull", "runif")
# distributionspparameters are the parameters for the specific function desired for speciation. 
# IMPORTANT: this vector of fuction parameters must *start by the second one, since the first parameter will always be one for all the function and is added already by this function*. HINT: see the help of the desired function for more details (e.g. ?rexp) Example of parameter for a exponential distribution with lambda of one (distributionspparameters <- c(1)). Entry in the distributionparameters can be "#", # or c(#,#) in case of more characters
# distributionextname is the same as the distributionspname but for the probability of extinction (e.g. distributionextname <- "rexp") 
# distributionextparameters is the same as the distributionspparameters but for the extinction probability function. By default extinction is set to ZERO, i.e. no extinction (e.g. distributionextparameters <- c(0)). Entry in can be "#", # or c(#,#) in case of more characters
# symmetric tells which macroevolutionary model should be used. If symmetric=TRUE the symmetric model will be used, else if FALSE, asymmetric model will be used. By default symmetric=TRUE
# complete: If complete = TRUE, the tree with the extinct and non-sampled lineages is returned. If complete = FALSE, the extinct and non-sampled lineages are suppressed. Complete=FALSE by default
# labellivingsp is the label that will be drawn on each tip surviving until the present. An automatic sequential number will be added to the chosen name. By default labellivingsp="sp."
# labelextinctsp is the label that will be drawn on each extinct tip. By default labelextinctsp <- "ext."
{
if (symmetric == TRUE) 
{
	mytree <- lapply(rep(age,numbsim), mytree.symmetric.age, distributionspname, distributionspparameters, distributionextname, distributionextparameters, complete, labellivingsp, labelextinctsp )
}
else
{
	mytree <- lapply(rep(age,numbsim), mytree.asymmetric.age, distributionspname, distributionspparameters, distributionextname, distributionextparameters, complete, labellivingsp, labelextinctsp )
}
}
return(mytree)	
}
