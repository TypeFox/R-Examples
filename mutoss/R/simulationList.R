# Generates a list of lists. Every of these lists, denoted for now
# by L, can be evaluated as a function call by eval(as.call(L))
# listFunAndParameter MUST have the following form:
# list(
#		funName = "UnifRandom",					# Description/Label of the function to be used
#		fun		= runif,					# A real function, not only the name
#		n		= 1						# "n" is integer   
# )
# Example for the argument:
# Arguments of runif are: n, min, max.
# listFunAndParameter = list(funName="UnifRandomVariable", fun=runif, n=2, min=c(1:2), max=c(1.1, 2.1))
# sapply(generateFunctionStack(listFunAndParameter), function(fc) eval(as.call(fc)))
# is same as
# runif(2, 1, 1.1); runif(2, 2, 1.1); runif(2, 1, 2.1); runif(2, 2, 2.1)
# Of course the second and third call do not make sense.
generateFunctionStack <- function(listFunAndParameter) 
{
	# special case: the procedure does not need any further parameter
	if (length(listFunAndParameter) == 2) 
		return(list(noParameters=list(fun=listFunAndParameter$fun)))
	
	# listFunAndParameter[[3]] is the first real parameter. The first 2 are the function to
	# to be called and a description.
	outerPar <- 1:length(listFunAndParameter[[3]])
	
	# Actually I want to build the outerproduct of the parameters,
	# but instead of this I use index numbers indicating the position
	# of the used parameter. If  n=c("a", "b"), n0=c(1:5), alpha=c(0.1, 0.2)
	# then 2, 3, 2 stands for n="b", n0=3, alpha=0.2
	for(par in listFunAndParameter[-3:-1])
		outerPar <- outer(outerPar, 1:length(par), paste)
	
	# special case of only ONE parameter
	if (length(listFunAndParameter) == 3)
		outerPar <- outer(outerPar, "", paste)		
	
	# build now for every parameter constellation
	# a list that can be casted into a function call. 
	fcStack <- list()
	for (parIDX in outerPar)
	{				
		idx <- as.numeric(unlist(strsplit(parIDX, " ")))
		
		parameter <- list()
		for(i in 1:length(listFunAndParameter[-2:-1]))
			# listFunAndParameter[-2:-1][[i]] is the i-th parameter in the list.
			# from the 1st parameter we want the idx[1]-th entry from the 2nd parameter
			# we want the idx[2]-th entry and so on.
			parameter <- c(parameter, listFunAndParameter[-2:-1][[i]][idx[i]])
		
		stackPosName <- paste(listFunAndParameter$funName, parIDX)
		fcStack[[stackPosName]] <- c(listFunAndParameter$fun, parameter)
		names(fcStack[[stackPosName]]) <- names(listFunAndParameter[-1])
	}	
	return(fcStack)
}


gatherParameters <- function(simObject) 
{
	#+++++++++++++++++++	Subfunctions	+++++++++++++++++++++++
	# extract from resultVecotr ( = simObject$results ) all values of the
	# parameter with the name paraName.
	
	getParamWithName <- function(resultVector, paraName)
	{	
		unlist( 
				lapply(resultVector,
						function(mts)
						{
							val <- mts$parameters[[paraName]]
							if (is.null(val)) 
								return("")
							
							val				
						}
				)
		)
	}	
	#-------------------	Subfunctions	-----------------------
	
	# gathering all parameters used in the simObject$results
	parNames <- unique(unlist(lapply(simObject$results, function(obj) names(obj$parameters))))
	
	return(data.frame(sapply(parNames, function(pN) getParamWithName(simObject$results, pN))))
	
#	# calling data.frame(sapply(parNames, function(pN) getParamWithName(simObject$results, pN)))
#	# is not good, because numeric parameters can be converted to characters and then to factors
#	# and it is possible that the original order is lost. For example the order of c(64, 128)
#	# will be 128 < 64.

#	ret = data.frame(getParamWithName(simObject$results, parNames[1]))
#
#	for (pN in parNames[2:length(parNames)]) # {"funName", "method"} subset of parNames; that is length(parName) >= 2 
#		ret = data.frame(ret, factor(getParamWithName(simObject$results, pN)))
#	
#	names(ret) <- parNames	
#	return(ret)
}


gatherStatistics <- function(simObject, listOfStatisticFunctions, listOfAvgFunctions) 
{
	#+++++++++++++++++++++++++++	Subfunctions	++++++++++++++++++++++++++
	# calculates the intersection of all elements given in aList
	listIntersect <- function(aList)
	{
		nn <- length(aList)		
		if(nn == 1)
			return(aList[[1]])
		
		intersect(aList[[1]], listIntersect(aList[-1]))
	}
	# actually the whole work is done by this subfunction.
	gatherStatisticsOneAvgFun <- function(simObject, listOfStatisticFunctions, avgFun, avgFunName = deparse(substitute(avgFun)))
	{
		# extract the parameter constellations form the obejct returned by simulation()
		paraNameDF <- gatherParameters(simObject)
		unqParaNameDF <- unique(paraNameDF)
		rownames(unqParaNameDF) <- 1:length(rownames(unqParaNameDF))
		
		# this will be the data.frame containing the parameter constellation and the calculated (averaged) statistics
		statDF <- data.frame()
		data.set.numbers <- sapply(simObject$results, function(res) res$data.set.number)
		for (i in rownames(unqParaNameDF))
		{
			# search which objects in simObject$results belong to 
			# parameter configuration in unqParaNameDF[i, ]
			idxs <- listIntersect(
					lapply(names(paraNameDF), 
							function(pN) which(unqParaNameDF[i, pN] == paraNameDF[ , pN])
					)
			)
			
			# applying any given statistic to the objects with the same 
			# parameter constellation			
			if (missing(avgFun))
			{ # no avgFun, thus the resulting data.frame will have one row for every simObject$results
				tmp <- sapply(listOfStatisticFunctions,	
						function(fun) sapply(idxs, function(idx) fun(simObject$data[[data.set.numbers[idx]]], simObject$results[[idx]]))
				)
				statDF <- rbind(statDF, cbind(paraNameDF[idxs,], tmp))
			}else
			{ 	# avgFun supplied, thus the resulting data.frame will have only one row for
				# every parameter constellation

				statDF <- rbind(statDF, 
						sapply(listOfStatisticFunctions,
								function(fun) avgFun(sapply(idxs, function(idx)
														fun(simObject$data[[data.set.numbers[idx]]],
																simObject$results[[idx]])
											)
									)
						)
				)
			}
		}
		
		if (missing(avgFun))
		{
			# number the rows consecutively
			rownames(statDF) <- 1:length(rownames(statDF))
			
			return(
					list(
							statisticDF = statDF, 
							name.parameters = names(paraNameDF), 
							name.statistics = names(listOfStatisticFunctions), 
							name.avgFun = ""
					)
			)
		}
		
		# label the columns of the resulting data.frame
		names(statDF) <- paste(names(listOfStatisticFunctions), avgFunName, sep=".")
		statDF <- cbind(unqParaNameDF, statDF)		
		
		list(
				statisticDF = statDF, 
				name.parameters = names(paraNameDF), 
				name.statistics = paste(names(listOfStatisticFunctions), avgFunName, sep="."), 
				name.avgFun = avgFunName
		)
	}
	#---------------------------	Subfunctions	--------------------------	
	
	# if no average function is given 
	# the resulting data.frame will have one row 
	# for every object in simObject$results
	if (missing(listOfAvgFunctions))
		return(gatherStatisticsOneAvgFun(simObject, listOfStatisticFunctions))
	
	# the average function is a function, pass this directly to 
	# gatherStatisticsOneAvgFun
	if (is.function(listOfAvgFunctions))
	{		
		return(gatherStatisticsOneAvgFun(
						simObject, 
						listOfStatisticFunctions, 
						listOfAvgFunctions, 
						deparse(substitute(listOfAvgFunctions))
				)
		)
	}
	
	# call gatherStatisticsOneAvgFun for every function in 
	# listOfAvgFunctions 
	if (length(listOfAvgFunctions) > 0)
	{		
		if (sum(names(listOfAvgFunctions) != "") != length(listOfAvgFunctions))
			warning("The functions in listOfAvgFunctions should have a name!")
		
		tmp <- list()
		# cnt is needed to determine the name of "fun"
		cnt <- 0
		for (fun in listOfAvgFunctions)
		{
			cnt <- cnt + 1
			tmp[[cnt]] <- gatherStatisticsOneAvgFun(
					simObject,
					listOfStatisticFunctions,
					fun,
					names(listOfAvgFunctions)[cnt]
			)
		}
		
		
		# We have gathered many statistics, now join the information 
		ret <- tmp[[1]]
		ret$statisticDF <- ret$statisticDF[ret$name.parameters]
		
		
		for (i in seq(along.with = listOfAvgFunctions))
		{
			ret$statisticDF <- cbind(ret$statisticDF, tmp[[i]]$statisticDF[tmp[[i]]$name.statistics])
			ret$name.statistics <- c(ret$name.statistics, tmp[[i]]$name.statistics)
			ret$name.avgFun <- c(ret$name.avgFun, tmp[[i]]$name.avgFun)
		}
		
		ret$name.statistics <- unique(ret$name.statistics)
		ret$name.avgFun <- unique(ret$name.avgFun)
		
		return(ret)
	}	
}


simulation <- function(replications, DataGen, listOfProcedures, discardProcInput=FALSE) 
{	
	paraNameDataGen <- names(DataGen)
	if (length(paraNameDataGen) != length(unique(paraNameDataGen)))
	{
		cat("Parameter of data generating function:\n\t", paraNameDataGen, "\n")
		stop("Parameternames of the data generating function are not unique")
	}
	
	# check if parameter of the procedures are unique
	nameProblems <- FALSE
	for( i in seq(along.with = listOfProcedures) )
	{
		paraNameProc <- names(listOfProcedures[[i]])
		if (length(paraNameProc) != length(unique(paraNameProc)))
		{
			nameProblems <- TRUE
			cat("Parameters of procedure", listOfProcedures[[i]]$funName, "are not unique:\n\t", paraNameProc, "\n")
		}
	}
	if (nameProblems) stop("Parameters of some procedures are not unique.\n")
	
	# no intersection between parameters of the data generating function and the
	# procedures are allowed.
	nameProblems <- FALSE	
	for( i in seq(along.with = listOfProcedures) )
	{
		paraNameProc <- names(listOfProcedures[[i]])
		equalNames = sort(intersect(paraNameDataGen, paraNameProc))
		
		if (length(equalNames)!= 2 || !all(equalNames == c("fun", "funName")))
		{
			nameProblems <- TRUE			
			cat("Common names of the data generating function and multiple test procedure", listOfProcedures[[i]]$funName, "are:\n\t", equalNames, "\n")
		}		
	}
	if (nameProblems) stop("The only common name of data generating function and multiple test procedure should be 'fun' and 'funName'.\n")
	
	# TODO: MS print progress of the simulation on the console!
	
	# generating all data generating functions
	dataGenStack 	<- generateFunctionStack(DataGen)
	
	# a bunch of stacks full of procedures
	# for example for every method (bonferroni and holm) there
	# is a stack for bonferroni with the different parameter configurations
	# and a stack for holm with the different parameter configurations
	procedureStacks	<- lapply(listOfProcedures, function(procs) generateFunctionStack(procs))
	names(procedureStacks) 	<- sapply(listOfProcedures, function(procs) procs$funName)
	
	ret = list()
	
	# cnt is used as an identifier. So every list with
	# the same $data.set.number is based on the same generated data
	cnt <- 0
	for( dataGenCall in dataGenStack )
	{	
		# This is probably the right place for gridComputation
		
		# It calls ONE time
		# dataGenCall. Every procedure in procedureStacks is applied
		# to this one "dataSet". returns $data and $results
		genOneDataSetAndApplyProcedures <- function(dummy) 
		{
			# generating data
			data <- eval(as.call(dataGenCall))
			
			# cnt is a global variable that has to be increased                     
			# each time a new "dataSet" is generated.
			assign("cnt", get("cnt", envir=sys.frame(-2)) + 1, envir = sys.frame(-2))
			
			# every procedure will be applied to the generated Dataset and the results
			# will be stored in the following list
			procs.results <- vector("list", sum(sapply(procedureStacks, function(stack) length(stack))))
			procs.results.idx <- 0
			for (pS in seq(along.with=procedureStacks))
			{
				# procedureStacks consists of stacks, go through one by one
				# this means applying every procedure to the given dataset
				procStack <- procedureStacks[[pS]]
				for (proc in procStack)
				{
					procs.results.idx <- procs.results.idx + 1
					result <- list()
					
					# every dataset get a unique number
					result$data.set.number <- get("cnt", envir  = sys.frame(-2))
					
					# saving the parameter constellation of the used data generating function
					paramDataGen <- c(DataGen$funName, dataGenCall[-1])
					names(paramDataGen)[1]	<- "funName"

					# saving the parameter constallation of the used procedure
					paramProc <- names(procedureStacks)[pS]
					
#					if (length(proc) == 1) 
#						#proc uses only the parameter from the output of the data generating function
#						paramProc <- c(paramProc)#, dummy="")
#					else
					if(length(proc)>1)
						paramProc <- c(paramProc, proc[-1])
					
					names(paramProc)[1]	<- "method"
					
					# saving the parameter constallation of the data generating function and the procedure 	
					result$parameters <- c(paramDataGen, paramProc)
					

					# next step is to assign the inputdata for the procedure parameters that was generated by
					# the data generating function. But at first I check if this will overwrite
					# other parameters already specified for the procedure.
					inter <- intersect(names(data$procInput), names(proc))
					if (length(inter) != 0)
						warning("\n\n\tSome of the parameter of one procedure are already specified,\n\t",
								"and the data generating function now provides new values for these",
								"parameters :\n\n\t",
								"Affected procedure : ", listOfProcedures[[pS]]$funName, "\n\t",
								"Affected parameters: ", paste(inter, collapse=" "), "\n")
					
					# assign inputdata generated by the data generating function to the procedure parameters.
					for(paraInputName in names(data$procInput))
						proc[[paraInputName]] <- data$procInput[[paraInputName]]
					
					# calling the procedure
					procOutput 	<- eval(as.call(proc))
					
					# writing the output of the procedure into the result
					for(name in names(procOutput))
						result[[name]] <-  procOutput[[name]]
					
					# append the new result
					procs.results[[procs.results.idx]] <- result
				}
			}

			if (discardProcInput)
				data$procInput = NULL		
			
			# return the used dataset, the output of the procedure with the used
			# parameter constellation
			return(list(data=data, results=procs.results))
		}
		ret <- c(ret, lapply(1:replications, genOneDataSetAndApplyProcedures))
	}
	

	# I want to have $data for data and $results for the output of
	# the procedures and parameter constellations. 	
	only.results <- sapply(ret, function(obj) obj$results)	
	dim(only.results) <- NULL
	return(
			list(data=lapply(ret, function(obj) obj$data),
					results=only.results
			)
	)	
}

