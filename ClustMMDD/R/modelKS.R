#################################################
##  Class modelKS of parameters: modelKS.R       ##
#################################################

##	@author	Wilson Toussile

##	Help
#'	The class of parameters for a given number K of components and a selected subset S of variables.
#'	@title	ModelKS is the class of parameters of a model (K, S).
#'	@slot	N	The size of the associated data.
#'	@slot	P	The number of variables.
#'	@slot	N_levels	The vector of the numbers of levels: N_levels[j] is the number of levels of variable j.
#'	@slot	levels	The list levels of the variables.
#'	@slot	K	The number of components in the mixture.
#'	@slot	S	The subset of relevant variables for clustering purposes.
#'	@slot	dim	The number of free parameters in the model defined by (K, S) (the dimension of the model).
#'	@slot	mixingProportions	The vector of the mixing proportions.
#'	@slot	count	The list of the observed counts of the levels for the variables.
#'	@slot	frequencies	The list of the observed frequencies of the levels for the variables.
#'	@slot	proba	The list of the proba of the levels for the K components.
#'	@slot	logLik	The maximum log-likelihood.
#'	@slot	entropy	The hard entropy.
#'	@slot	membershipProba	The matrix of the membership proba.
#'	@slot	mapClassification	The vector of membership assignment by MAP.
#'	@references	 W. Toussile and E. Gassiat (2009)
#'	@author	Wilson Toussile
# #' @examples
# #'	data(par_ks)
# #'	par_ks	
#'	@exportClass	modelKS


setClass(Class = "modelKS",
  slots = c(
	N = "integer", #	Size of data
	P = "integer", # Number of looci 
	N_levels = "integer", # Number of alleles per loci
	levels = "list",
	K = "integer",#  Number of classes
	S = "logical", # Vector of indicators of selected loci                 
	dim = "integer",  # Dimension = number of free parameters
	mixingProportions = "numeric",  # Mixing proportions
	count = "list", # List of matrix of allelic frequencies per locus
	frequencies = "list",	#	List of observed frequencies
	proba = "list",
	logLik = "numeric", # Log-likelihood
	entropy = "numeric",	#	Soft entropy
	membershipProba = "matrix", # Membership probabilities. It is a matrix of numerics.
	mapClassification = "integer"	# MAP classification
  ),

  ##  Validity TODO
  validity = function(object)
	{
	  ##  TODO
	  if(object@P>0 && length(object@S) != object@P) stop("[modelKS:validation] Invalid length of S !")
	  if(object@P>0 && length(object@N_levels) != object@P) stop("[modelKS:validation] Invalid length of N_levels !")
	  if(object@P>0 && length(object@levels) != object@P) stop("[modelKS:validation] Invalid length of levels !")
	  if(length(object@count)>0 && object@P>0 && length(object@count) != object@P) stop("[modelKS:validation] Invalid length of count !")
	  if(object@P>0 && length(object@proba) != object@P) stop("[modelKS:validation] Invalid length of proba !")
	  if(object@K && length(object@mixingProportions) != object@K) stop("[modelKS:validation] Invalid length of mixingProportions !")
	  
	  if(object@K>1 & object@P>0)
	  {
		for(j in 1:object@P) 
		{
		  if(ncol(object@proba[[j]]) != object@K) stop("[modelKS:validation] Invalid ncol in some proba[[j]] !")
		}
	  }
	  
	  if(object@P>0 & object@K>0)
	  {
		if(object@K>1)
		{
		  for(j in 1:object@P)
		  {
			if(nrow(object@proba[[j]]) != object@N_levels[j]) stop("[modelKS:validation] Invalid ncol in some proba[[j]]")
		  }
		}
		else
		{
		  for(j in 1:object@P)
		  {
			if(length(object@proba[[j]]) != object@N_levels[j]) stop("[modelKS:validation] Invalid length of some proba[[j]]")
		  }
		}
	  }
	  cat("***  End modelKS:validity method\n")
	  return(TRUE)
	}
)

##  Show method	##
#'	Show an object of class [\code{\linkS4class{modelKS}}].
#'	@title	show
#'	@param	this	An object of class [\code{\linkS4class{modelKS}}].
#'	@author	Wilson Toussile
# #'	@exemples
# #'	data(par_ks)
# #'	par_ks	
#'	@export

setMethod(f="show", signature = "modelKS",
  function(object)
  {
    cat("** Print a set of paramters of modelKS class **\n")
	cat("\nSize of the dataset N = ", object@N)
    cat("\nNumber of variables P = ", object@P)
    cat("\nThe numbers of observed levels N_levels = ", object@N_levels)
	
	cat("\n\n** Model (K, S) : \n")
    cat("\nNumber of clusters K = ", object@K)
	cat("\nClustering variables S = ", object@S)

    cat("\n\n** Mixing proportions : \n")
    cat("Mixing proportions mixingProportions = ", object@mixingProportions) 
	
    cat("\n\nProbabilities in clusters : \n"); print(object@proba)
    
    cat("\nNumber of free parameters dim = ", object@dim)
    cat("\nLog likelihood = ", object@logLik)
    cat("\nEntropy = ", object@entropy)
	
	#cat("\n\nLevels count : \n"); print(object@count)
	
	#cat("\n\nObserved frequencies : \n"); print(object@frequencies)
	
    cat("\n\n***  End modelKS:show method ***\n")
  }
)

## Set from a list
setGeneric( name = "setModelKS", package = "ClustMMDD", 
  def = function(x)
  {
	standardGeneric("setModelKS")
  }
)

## Set from a list
setMethod(f = "setModelKS", signature = "list",
  definition = function(x = "list")
  {	
	
	S = x[["S"]]
	K = as.integer(x[["K"]])
	for(i in 1:length(x[["prob"]]))
	{
	  rownames(x[["prob"]][[i]]) = x[["levels"]][[i]]
	  colnames(x[["prob"]][[i]]) = paste("Cluster", 1:K, sep="")
	}
	
	new(Class = "modelKS", N = as.integer(x[["N"]]), 
		P = as.integer(length(S)), 
		N_levels = as.integer(x[["NbersLevels"]]), 
		levels = x[["levels"]], 
		K = K, 
		S = S, 
		dim = as.integer(x[["dim"]]), 
		mixingProportions = as.double(x[["pi_K"]]), 
		proba = x[["prob"]], 
		logLik = as.double(x[["logLik"]]), 
		entropy = as.double(x[["entropy"]]), 
		membershipProba = x[["Tik"]], 
		mapClassification = as.integer(x[["mapClassif"]]))
  }
)


##	Replace a slot ##
#'	Method that alows replacing a value of a slot for class [\code{\linkS4class{modelKS}}].
#'	@title	Operator []
#'	@param	...	The name of order of the slot to be replace.
#'	@author	Wilson Toussile
#'	@export
setMethod(f = "[<-", signature = "modelKS",
	definition = function(x, i, j, value)
	{
		switch(EXPR=i,
			"N" = {x@N<-value[1]},
			"P" = {x@P<-value[1]},
			"N_levels" = {x@N_levels<-value},
			"levels" = {x@levels<-value},
			"K" = {x@K<-value},
			"S" = {x@S<-value},
			"mixingProportions" = {x@mixingProportions<-value},
			"dim" = {x@dim<-value},
			"count" = {x@count<-value},
			"frequencies" = {x@frequencies<-value},
			"proba" = {x@proba<-value},
			"logLik" = {x@logLik<-value},
			"entropy" = {x@entropy<-value},
			"membershipProba" = {x@membershipProba<-value},
			"mapClassification" = {x@mapClassification},
			stop("This attribute does not exist !")
		)
		return(x)
	}
)


##	Get a slot ##
#'	Method that alows replacing a value of a slot for class [\code{\linkS4class{modelKS}}].
#'	@title	operator []
#'	@param	...	The name of order of the slot to be replace.
#'	@author	Wilson Toussile
# #'	@examples
# #'	data(par_ks)
# #'	par_ks	
# #		par_ks["PROB"]
#'	@export

setMethod(f = "[", signature = "modelKS",
	definition = function(x, i, j, value)
	{
		switch(EXPR=i,
		"N" = {return(x@N)},
		"P" = {return(x@P)},
		"N_levels" = {return(x@N_levels)},
		"levels" = {return(x@levels)},
		"K" = {return(x@K)},
		"S" = {return(x@S)},
		"mixingProportions" = {return(x@mixingProportions)},
		"dim" = {return(x@dim)},
		"count" = {return(x@count)},
		"frequencies" = {return(x@frequencies)},
		"proba" = {return(x@proba)},
		"logLik" = {return(x@logLik)},
		"entropy" = {return(x@entropy)},
		"membershipProba" = {return(x@membershipProba)},
		"mapClassification" = {return(x@mapClassification)},
		stop("This attribute does not exist !")
		)
		return(x)
	}
)


##	Operator == 
#'	Compare two models by returning the indicator of if the have the same number of components and the same subset of relevant variavles.
#'	@title	==
#'	@param	e1	An object of class [\code{\linkS4class{modelKS}}].
#'	@param	e2	An object of class [\code{\linkS4class{modelKS}}].
#'	@return	TRUE if K(e1) == K(e2) and S(e1) == S(e2), and FALSE otherwise.
#'	@author	Wilson Toussile
#'	@export

setMethod(f = "==", signature = "modelKS",
  definition = function(e1, e2)
  {
	#	par1 and par2 are of type modelKS
	S1 = e1["S"]
	S2 = e2["S"]
	K1 = e1["K"]
	K2 = e2["K"]
	if(((K1==1) | (all(!S1))) && ((K2==1) | (all(!S2)))) return(TRUE)
	return(length(S1) == length(S2) && K1 == K2 && all(S1 == S2))
  }
)

##	is.element ##
#'	Return TRUE if an object of class [\code{\linkS4class{modelKS}}] belong to a set of objects of the same class.
#'	@title	is.element
#'	@param	el	An object of class [\code{\linkS4class{modelKS}}].
#'	@param	set	A set of objects of class [\code{\linkS4class{modelKS}}].
#'	@return	TRUE if the object el belongs set.
#'	@author	Wilson Toussile
#'	export


setMethod(f = "is.element", signature(el="modelKS", set="modelKS"),
  definition = function(el = "modelKS", set = "modelKS")
  {	
	el = c(el)
	set = c(set)
	return(any(sapply(set, function(xpar) {return(el == xpar)})))
  }
)


#########################################
##	Read modelKS from a file
setGeneric( name = "read.modelKS", package = "ClustMMDD", 
  def = function(file)
  {
	standardGeneric("read.modelKS")
  }
)

setMethod(f = "read.modelKS", signature("character"),
  definition = function(file)
  {
	xmodelKS = readParKS_Rcpp(file)
	xmodelKS[["Tik"]] = matrix() # Set membershipProba to an empty matrix
	return(setModelKS(xmodelKS))
  }
)

#########################################
##	Read modelKS from a file
setGeneric( name = "model", package = "ClustMMDD", 
  def = function(object)
  {
	standardGeneric("model")
  }
)

setMethod(f = "model", signature("modelKS"),
  definition = function(object)
  {
	return(list(K = object@K, S = object@S))
  }
)


#########################################
##	Write in a file
# setGeneric( name = "write.modelKS", package = "ClustMMDD", 
#   def = function(object, file)
#   {
# 	standardGeneric("write.modelKS")
#   }
# )
# 
# 
# 
# setMethod(f = "write.modelKS", signature(object = "modelKS", file="character"),
#   definition = function(object = "modelKS", file = "character")
#   {	
# 	## TODO
# 	###writeParKS_InFile_Rcpp(Rcpp::List modelList, std::string file)
#   }
# )


#########################################
##	Is an object is from modelKS?

setGeneric( name = "is.modelKS", package = "ClustMMDD", 
  def = function(object)
  {
	standardGeneric("is.modelKS")
  }
)

setMethod(f = "is.modelKS", signature("ANY"),
  definition = function(object)
  {
	return(is(object, "modelKS"))
  }
)


#########################################
##	Simulation method
#'	Simulate a dataset using the parameters in a [\code{\linkS4class{PAR_KS}}] object.
#'	@title	The method simul.dist.data of [\code{\linkS4class{PAR_KS}}] class.
#'	@param	this	An object of class [\code{\linkS4class{PAR_KS}}].
#'	@param	N	The size of the dataset to simulate.
#'	@param	ploidy	The number of occurrences for an observation of a variable. For example, ploidy = 2 for genotypic data from diploid individual.
#'	@return	A list: the first element is the dataset simulated, and the second the membership.
#'	@author	Wilson Toussile
#'	@export

setGeneric( name = "simulData", package = "ClustMMDD",
  def = function(object, N, ploidy)
  {
	standardGeneric("simulData")
  }
)



setMethod( f = "simulData", signature("modelKS", "numeric", "numeric"),
  definition = function(object, N, ploidy)
  {
	Z = sample(1:object@K, N, replace = TRUE, prob=object@mixingProportions)
	
	MaxPloidy = 5
	
	if(!is.element(ploidy, c(1:MaxPloidy)))
	  stop("Incorrect argument ploidy !")
	
	if(ploidy == 1)
	{
	  data_simulated = sapply(1:object@P, function(j)
	  {
		if(object@S[j])
		{
		  Vj = rep(0, N)
		  for(k in 1:object@K)
		  {
			v = (Z==k)
			Vj[v] = sample(object@N_levels[j], sum(v), replace=TRUE, prob=object@proba[[j]][, k])
		  }
		  return(Vj)
		}
		else return(sample(object@N_levels[j], N, replace=TRUE, prob=object@proba[[j]][,1]))
	  })
	}
	else
	{
	  data_simulated = sapply(1:object@P, function(j)
	  {
		V = matrix(0, N, ploidy)
		if(object@S[j])
		{
		  for(jj in 1:ploidy)
		  {
			for(k in 1:object@K)
			{
			  Zk = (Z==k)
			  V[Zk, jj] = sample(object@N_levels[j], sum(Zk), replace=TRUE, prob=object@proba[[j]][, k])
			}
		  }
		  return(V)
		}
		else 
		{
		  for(jj in 1:ploidy)
		  {
			for(k in 1:object@K)
			{
			  Zk = (Z==k)
			  V[Zk, jj] = sample(object@N_levels[j], sum(Zk), replace=TRUE, prob=object@proba[[j]][, 1])
			}
		  }
		  return(V)
		}
	  })
	  data_simulated = matrix(data_simulated, ncol=ploidy*object@P)
	}
  
	return(list(data = data_simulated, class = Z))
  })




