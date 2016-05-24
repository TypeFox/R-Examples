# <OWNER> = Saip CISS
# <ORGANIZATION> = QueensBridge Quantitative
# <YEAR> = 2014

# LICENSE 
# BSD 3-CLAUSE LICENSE

# Copyright (c) 2014, Saip CISS
# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# END OF LICENSE 

# remove characters
rm.string = function(string, which.char, exact = TRUE)
{
	string2 = vector()
	n = length(string)
	
	for (i in 1:n)
	{	
		string2[i] = strsplit(string[i], which.char, fixed = exact)	 
	
		if (length(string2[[i]]) == 2)
		{ string2[i] = paste(string2[[i]][1], string2[[i]][2], sep ="" ) }
	}
	
	return(unlist(string2))	
}

#concat
concatCore <- function(X)
{
	concatX = NULL
	for (i in 1:length(X))
	{	concatX =  paste(concatX, X[i], sep="") }
	return(concatX)
}

concat <- function(X) apply(X, 1, function(Z) concatCore(as.character(Z)))

# most frequent values
modX <- function(X) as.numeric(names(which.max(table(round(X,0)))))

# indexes
find.idx =  function(sample_ID_vector, ID_vect, sorting = TRUE, replacement = 0)
{
	if (sorting == TRUE)
	{	ghost_idx_file = match(sort(ID_vect), sort(sample_ID_vector))	}
	else
	{	ghost_idx_file = match(ID_vect, sample_ID_vector)	}
	
	return(which(!is.na(ghost_idx_file)))	
}

find.first.idx <-  function(X, sorting = TRUE)
{
	all_idx = which.is.duplicate(X, sorting = sorting)
	first_idx = unique(match(X[all_idx],X)) 
	return(first_idx)
}

which.is.duplicate <- function(X, sorting = TRUE) 
{
	if (sorting)  {	 X =sort(X)	}
	return(na.omit(which(duplicated(X) == TRUE)) )
}

# FACTORS AND CHARACTERS
factor2vector = function(d)
{
	if (is.factor(d))
	{
		Lev = levels(d);
		nb_Lev = length(Lev)
		{ 
			G = vector(length = length(d));
			if (is.numeric(Lev[1]))
			{
				if (as.numeric(Lev[1]) == 0)
				{	G = d-1	}
			}
			else
			{ 	G =cbind(d)	}
			d = G;
		}
		return(list(vector = as.numeric(d), factors = Lev) )
	}
	else
	{ return(list(vector = as.vector(d), factors = NA)) }
}

which.is.factor <- function(X, maxClasses = floor(0.01*min(3000, nrow(X))+2), count = FALSE) 
{
	p = ncol(X)
	n = nrow(X)
	values = rep(0,p)
	options(warn = -1)
	
	for (j in 1:p) 
	{ 
		flag = FALSE
		uniqueValues = unique(X[,j])
		countValues <- length(uniqueValues)
				
		if (is.factor(X[,j]))
		{ 
			XLevels = levels(uniqueValues)
			checkThis = which(XLevels == "?")
			if (length(checkThis) > 0) { XLevels =  XLevels[-checkThis] }

			checkNA = which(is.na(XLevels))
			if (length(checkNA) > 0) { XLevels =  XLevels[-checkNA] }
			
			factorAsNumeric <- rmNA(as.numeric(XLevels))
			
			if (length(factorAsNumeric) == 0) { flag = TRUE }
			else
			{
				if (!is.double(factorAsNumeric))
				{ 
					if ( (max(factorAsNumeric) - min(factorAsNumeric)) == (countValues - 1) )	{	flag  = TRUE  }
					
					if (countValues <= maxClasses) 	{  flag = TRUE }
				}
			}			
						
			if (flag) 
			{
				if (count) { values[j] = countValues }
				else {  values[j] = 1 	}	
			}
		}
		else
		{
			if (!is.numeric(X[,j]))
			{  
				if (countValues <= maxClasses)
				{ 
					if (count) { values[j] = countValues  }
					else { 	values[j] = 1 }	
				}
			}
		} 	
	}
	names(values) = colnames(X)
	options(warn = 0)
	
	return(values)
}


# MATRIX
factor2matrix = function(X, threads = "auto")
{
	if (is.factor(X) | is.vector(X))
	{ return(as.true.matrix(factor2vector(X)$vector)) }
	else
	{
		np = dim(X)
		if  (max(np) < 1000)
		{			
			for (j in 1:np[2])
			{	
				if (is.character(X[,j]))	{ X[,j] = as.factor(X[,j]) }
				
				X[,j] = factor2vector(X[,j])$vector	
			}
			return(as.true.matrix(X))
		}
		else
		{
			#require(parallel)	
			max_threads = detectCores(logical = FALSE)	
			
			if (threads == "auto")
			{	
				if (max_threads == 2) { threads = max_threads }
				else {	threads  = max(1, max_threads - 1)  }
			}
			else
			{
				if (max_threads < threads) 
				{	
					maxLogicalThreads = detectCores()
					if (maxLogicalThreads < threads)
					{ cat("Warning : number of threads indicated by user was higher than logical threads in this computer.\n") }
				}
			}
			
			#require(doParallel)				
			Cl = makePSOCKcluster(threads, type = "SOCK")
			registerDoParallel(Cl)
			chunkSize  <-  ceiling(np[2]/getDoParWorkers())
			smpopts  <- list(chunkSize  =  chunkSize)
			Z = matrix(NA, np[1], np[2])			
			Z <- foreach(j =1:np[2], .export = "factor2vector", .options.smp = smpopts, .combine = cbind) %dopar%
			{
				if (is.character(X[,j]))	{ X[,j] = as.factor(X[,j]) }
				factor2vector(X[,j])$vector	
			}
			
			colnames(Z) <- colnames(X)
			stopCluster(Cl)
			
			return(as.true.matrix(Z))
		}
	}
}

matrix2factor2 <- function(X, maxClasses = 10)
{
	idx = which.is.factor(X,maxClasses = maxClasses)
	factors = which(idx != 0)
	if (length(factors > 0))
	{	for (i in 1:length(factors)) {	X[,factors[i]] = as.factor(X[,factors[i]])	}	}
	else 
	{ cat("No variable found as categorical. Please try to increase number of classes.\n") }
	return(X)
}


NAfactor2matrix <- function(X, toGrep = "")
{
	if (length(toGrep) > 1)
	{
		for (i in 1:length(toGrep))
		{  
			toGrepIdx = NULL
			toGrepIdx = rbind(toGrepIdx, which(X == toGrep[i], arr.ind = TRUE))				
		}
	}
	else
	{	
		if (toGrep == "anythingParticular")	{ toGrepIdx = integer(0) }
		else  { toGrepIdx = which(X == toGrep, arr.ind = TRUE)  }		
	}
	
	NAIdx = which(is.na(X), arr.ind = TRUE)
	
	X <- factor2matrix(X)
	
	if (length(dim(toGrepIdx)[1]) > 0) { X[toGrepIdx] = NA	}
	if (length(dim(NAIdx)[1]) > 0) {  X[NAIdx] = NA	  }
	
	return(X)
}

NAFeatures <- function(X) apply(X,2, function(Z) { if (length(rmNA(Z)) == 0) { 1 } else{ 0 } })

rmNA = function(X)
{
	NAIdx = which(is.na(X))	
	if (length(NAIdx) > 0) { return(X[-NAIdx])	}
	else { return(X) }
}

NATreatment <- function(X, Y, na.action = c("fastImpute", "accurateImpute", "omit"), regression = TRUE)
{
	p = ncol(X)
	featuresWithNAOnly = which(NAFeatures(X) == 1)
	if (length(featuresWithNAOnly) == p) { stop("Data have only NA values") }
	if (length(featuresWithNAOnly) > 0) 
	{ 
		X = X[,-featuresWithNAOnly] 
		cat("\nVariables with NA values only have been found and removed.\n")
	}
	NAvalues = which(is.na(X), arr.ind = TRUE)
	NAInputs = NAvalues[,1]
	matchNA = (length(NAInputs) > 0)		
	if ( (length(unique(NAInputs)) > (nrow(X) - 30)) & (na.action[1] == "omit") ) 
	{ stop("Too much missing values in data. Please impute missing values in order to learn the data.\n") }		
	if (is.factor(Y)) levelsY = levels(Y)
	if (!regression & !is.factor(Y) ) { Y = as.factor(Y); levelsY = levels(Y) }	
	
	NALabels = which(is.na(Y))
	matchNALabels = (length(NALabels) > 0)		
	if ( (length(NALabels) > (length(Y) - 30)) & (na.action[1] == "omit") )	{ stop("Too much missing values in responses.\n") }			
	if (matchNA | matchNALabels)
	{
		if (!is.null(Y))
		{
			if (matchNALabels) { newY = Y }
			else  { newY <- as.vector(NAfactor2matrix(as.numeric(Y))) }
			if (na.action[1] != "omit") 
			{ 
				cat("NA found in data. Imputation (fast or accurate, depending on option) is used for missing values.\nIt is strongly recommended to impute values outside modelling, using one of many available models.\n")					
				XY <- na.impute(cbind(X, newY), type = na.action[1])
				Y <- XY[,p+1]
				X <- XY[,-(p+1)]
				if (!regression) {  Y = as.factor(Y); levels(Y) = levelsY }
				rm(XY)				
			}
			else
			{
				if (matchNALabels & matchNA) {	rmRows = unique(c(NAInputs, NALabels))	}
				else
				{
					rmRows = NULL
					if (matchNA) {	rmRows = c(rmRows, unique(NAInputs))	}						
					if (matchNALabels) 	{	rmRows = c(rmRows,unique(NALabels))	}
				}					 
				X = X[-rmRows,]
				Y = Y[-rmRows]
			}
		}
		else
		{
			cat("\nIf accuracy is needed, it is strongly recommended to impute values outside modelling, using one of many available models.\n")
			if (na.action[1] != "omit") { 	X <- na.impute(X, type = na.action[1])	}
			else
			{		
				if (matchNA) { X <- X[-NAInputs,] 	}		
			}
		}
	}
		
	return(list(X = X, Y = Y))
}		


rmInf = function(X)
{
	InfIdx = which((X == Inf) | (X == -Inf))	
	if (length(InfIdx) > 0) { return(X[-InfIdx])	}
	else { return(X) }
}	

vector2factor = function(V) as.factor(V)


matrix2factor = function(X, which.columns) 
{
	X =	as.data.frame(X)
	for (j in which.columns)
	{	X[,j] = as.factor(X[,j])	}
	
	return(X)
}

vector2matrix <- function(X, nbcol) if (is.vector(X)) { matrix(data = X, ncol = nbcol) } else { X }

fillVariablesNames <- function(X)
{
	if(is.null(colnames(X))) 
	{  
		varNames = NULL
		for (i in 1:ncol(X)) { varNames = c(varNames,paste("V", i, sep="")) }
		colnames(X) = varNames
	}
	
	return(X)
}

as.true.matrix = function(X)
{
	if (is.matrix(X))
	{ return(X) }
	else
	{
		if (is.vector(X))
		{	
			names_X = names(X)
			X = t(as.numeric(as.vector(X)))
			if (!is.null( names_X ) ) { colnames(X) = names_X	}		
		}
		else
		{
			if (is.data.frame(X))
			{	
				names_X = colnames(X)
				X = as.matrix(X)  	
				col_X = ncol(X)				
				X = mapply( function(i) as.numeric(factor2vector(X[,i])$vector), 1:col_X)
				if (is.vector(X)) { X = t(X) }				
				if (!is.null(names_X)) { colnames(X) = names_X	}		
			}		
		}
		return(X)
	}
}

sortMatrix = function(X, which.col, decrease = FALSE)
{
	if (is.character(which.col))
	{ 	which.col = which( colnames(X) == which.col)	}
	
	return(X[order(X[,which.col], decreasing = decrease),])	
}

sortDataframe <- function(X, which.col, decrease = FALSE)
{
	if (is.character(which.col))
	{ 	which.col = which( colnames(X) == which.col)	}

	idx = order(X[,which.col], decreasing = decrease)
	
	return(X[idx,])
}	
	
genericCbind <- function(X,Y, ID = 1)
{
	matchIdx = match(X[,ID],Y[,ID])
	
	if (all(is.na(matchIdx)))
	{  	stop("Matrix do not share same values in 'ID'.\n") 	}
	else
	{
		naIdx = which(is.na(matchIdx))
		
		Z = matrix(data = 0, ncol = ncol(Y), nrow = nrow(X))
		colnames(Z) = colnames(Y)
		
		Z[which(!is.na(matchIdx)),] = Y[na.omit(matchIdx),]
		Z[naIdx,ID] = X[naIdx,ID]
		
		newZ = cbind(X,Z[,-1])

		return(newZ)
	}
}
	
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

which.is.wholenumber <- function(X)  which(is.wholenumber(X))

fillWith <- function(X, toPut = "NA")
{
	for (j in 1:ncol(X))
	{ 
		if ( is.factor(X[,j]) | is.character(X[,j]) )
		{
			X[,j]= as.character(X[,j])
			otheridx = which(X[,j] == "")
			
			if (length(otheridx) > 0)
			{
				levelsLength = length(levels(X[,j]))
				levels(X[,j])[levelsLength]= toPut 
				X[otheridx,j] = toPut
			}
		}
	}
	
	return(X)
}

na.impute <- function(X, type = c("fastImpute", "accurateImpute"))
{
	if (type[1] == "fastImpute")
	{	return(na.replace(X, fast = TRUE))	}
	else
	{ 	return(fillNA2.randomUniformForest(X)) 	}
}


# sample replacement (good and automatic, but prefere na.impute for accuracy)
na.replace = function(X, fast = FALSE)
{
	na.matrix = which(is.na(X), arr.ind = TRUE)
	n = nrow(X)
	factors <- which.is.factor(X)
	
	if (is.null(dim(na.matrix))) 
	{ return(X)	}
	else
	{ 
		p = unique(na.matrix[,2])
		for (j in seq_along(p))
		{
			i.idx = na.matrix[which(na.matrix[,2] == p[j]),1]
			
			if (!fast)
			{
				subLength <- length(rmNA(X[,p[j]]))
				for (i in seq_along(i.idx))
				{  
					if (factors[p[j]] == 1)	{  X[i.idx[i],p[j]] <- modX(sample(rmNA(X[,p[j]]),subLength, replace = TRUE)) 	}
					else
					{
						while(is.na(X[i.idx[i],p[j]])) 	{ 	X[i.idx[i],p[j]] <-  mean(rmNA(X[sample(1:n,sample(1:n,n)),p[j]])) 	}
						
						if ( length(which.is.wholenumber(rmNA(X[,p[j]]))) == subLength ) 	
						{  X[i.idx[i],p[j]] =  floor(X[i.idx[i],p[j]]) }
					}
				}
			}
			else
			{	
				if (factors[p[j]] == 1)
				{ 	X[i.idx ,p[j]] = modX(rmNA(X[,p[j]]))	}
				else
				{	X[i.idx,p[j]] = median(rmNA(X[,p[j]]))	}
			}
		}
		
		return(X)
	}
}


parallelNA.replace <- function(X, threads = "auto", parallelPackage = "doParallel", fast = FALSE)
{
	p = ncol(X)
	# code parallelization
	{
		export = "na.replace"	
		#require(parallel)		
		max_threads = detectCores(logical = FALSE)		
		if (threads == "auto")
		{	
			if (max_threads == 2) { threads = max_threads }
			else {	threads = max(1, max_threads - 1)  }
		}
		else
		{
			if (max_threads < threads) 
			{	
				maxLogicalThreads = detectCores()
				if (maxLogicalThreads < threads)
				{ 
					cat("Warning : number of threads indicated by user was higher than logical threads in this computer.") 
				}
			}
		}
				
		{
			#require(doParallel)			
			Cl = makePSOCKcluster(threads, type = "SOCK")
			registerDoParallel(Cl)
		}
	}
		
	X <- foreach(i = 1:p, .export = export, .inorder = TRUE, .combine = cbind) %dopar%
	{  na.replace(X, fast = fast)	}
	
	return(X)
}


pseudoNAReplace = function(X, toGrep = "NA")
{
	
	p =ncol(X)
	for (j in 1:p)
	{
		if( ( is.character(X[,j]) | is.factor(X[,j]) ) )
		{
			X[,j] = as.character(X[,j])
			idx = which(X[,j] == toGrep)
			
			if (length(idx) > 0)
			{
				factorTable = table(X[,j])
				factorTableNames = names(factorTable)
				factorTableNames = factorTableNames[-which(factorTableNames == toGrep)]
				
				sampleIdx = sample(seq_along(factorTableNames),length(idx), replace = TRUE)
				X[idx, j] = factorTableNames[sampleIdx]
			}
		}
	}
	
	return(X)
}
	

na.missing = function(X, missing.replace = 0)
{
	na.data = which(is.na(X))
	X[na.data] = missing.replace
	return(X)
}

which.is.na = function(X) 
{
	if (!is.vector(X))
	{ 	X = factor2vector(X)$vector }
	
	return(which(is.na(X)))
}

# random combination : only multiple of 2 variables
randomCombination <- function(X, combination = 2, weights = NULL)
{
	howMany = ifelse(length(combination) > 1, length(combination), combination)
	n = nrow(X); p = ncol(X)
	L.combination = length(combination)
	if (L.combination > 1)
	{  
		var.samples = combination
		if (is.null(weights))
		{	var.weights = replicate(L.combination/2, sample(c(-1,1),1)*runif(1))	}
		else
		{	var.weights =  weights }
			
		X.tmp = matrix(NA, nrow(X), L.combination/2)
		colnames(X.tmp)= rep("V", ncol(X.tmp))
		idx = 1
		for (j in 1:(L.combination/2))
		{	
			X.tmp[,j] = var.weights[j]*X[,combination[idx]] + (1 - var.weights[j])*X[,combination[idx+1]]
			idx = idx + 2
		}
		
		idx = 1
		for (i in 1:ncol(X.tmp)) 
		{ 
			colnames(X.tmp)[i] = paste("V", combination[idx], "x", combination[idx+1], sep="") 
			idx = idx + 2
		}
		X.tmp = cbind(X, X.tmp)
		
	}
	else
	{ 
		var.samples = sample(1:p, 2*howMany)	
		var.weights = replicate(howMany, sample(c(-1,1),1)*runif(1))
		var.weights = cbind(var.weights, 1 - var.weights)

		X.tmp = X
		idx = 1
		for (j in 1:(length(var.samples)/2))
		{	
			X.tmp[,var.samples[idx]] = var.weights[j,1]*X[,var.samples[j]] + var.weights[j,2]*X[,var.samples[idx+1]]
			idx = idx + 2
		}
	}
	
	return(X.tmp)	
}
	

# LISTS
mergeLists <- function(X, Y, add = FALSE, OOBEval = FALSE)
{
	n = length(X)
	if (length(Y) != n)
	{ stop("lists size not equal.\n") }
	else
	{ 	
		Z = vector('list', n)
		names(Z) = names(X)
		for (i in 1:n)
	    {    
			if (add)
			{ 	Z[[i]] = X[[i]] + Y[[i]]  }
			else
			{   
				if (OOBEval)
				{
					if (is.matrix(X[[i]]))
					{	
						pX = ncol(X[[i]]); pY = ncol(Y[[i]])
						if ( pX == pY  )	
						{	
							if (nrow(X[[i]]) == nrow(Y[[i]])) {	Z[[i]] = cbind(X[[i]],Y[[i]])	}
							else {	Z[[i]] = rbind(X[[i]],Y[[i]])	}
						}
						else
						{
							if ( nrow(X[[i]]) == nrow(Y[[i]]) )	{	Z[[i]] = cbind(X[[i]],Y[[i]])	}
							else 
							{	
								pSamples <- (pX - pY)
								
								if (pSamples < 0)
								{	
									newM <- matrix(apply(X[[i]], 1, function(Z) sample(Z, abs(pSamples), replace= TRUE)), ncol = abs(pSamples))
									XX = cbind(X[[i]],newM)
									YY = Y[[i]]
								}
								else
								{
									newM <- matrix(apply(Y[[i]], 1, function(Z) sample(Z, pSamples, replace= TRUE)), ncol = pSamples)
									YY = cbind(Y[[i]], newM)
									XX = X[[i]]
								}								
															
								ZZ = rbind(XX,YY)
								newM2 <- matrix(apply(ZZ, 1, function(Z) sample(Z, min(pX, pY), replace= TRUE)), ncol = min(pX, pY))
							
								Z[[i]] = cbind(ZZ, newM2)
							
							}
						}							
					}
					else
					{	
						if ( (is.null(X[[i]])) & (is.null(Y[[i]])) ) { Z[[i]] = NULL }
						else {	Z[[i]] = c(X[[i]],Y[[i]])	}
					}
				}
				else
				{
					if (is.matrix(X[[i]])) 	{	Z[[i]] = rbind(X[[i]],Y[[i]])  }
					else 
					{		
						if ( (is.null(X[[i]])) & (is.null(Y[[i]])) ) { Z[[i]] = NULL }
						else {	Z[[i]] = c(X[[i]],Y[[i]])	}					
					}
				}
			}
		}
		names(Z) = names(X)
		return(Z)
	}
}


rm.InAList = function(X,idx)
{ 	
	tmp.X = X
	nb.list =length(X)
	new.X = vector('list',nb.list - length(idx))
	j = 1
	
	for (i in (1:nb.list)[-idx])
	{ 
		new.X[[j]] = X[[i]]
		names(new.X)[j] = names(X)[i] 
		j = j + 1
	}
	
	return(new.X)
}

rmInAListByNames <- function(X, objectListNames)
{
	for (i in 1:length(objectListNames))
	{
		idx = which(names(X) == objectListNames[i])
		
		if (length(idx) > 0)  { X = rm.InAList(X,idx) }
	}
	
	return(X)
}
		

# INSERT
insert.in.vector = function(X, new.X, idx)
{
	Y = vector(length = (length(X) + length(idx)) )
	n = length(Y)
	flag = 0;
	if ( idx <= (length(X) + 1) )
	{
		if (idx == (length(X) + 1))
		{	Y = c(X,new.X)	}
		else
		{		
			for (i in 1:n)
			{
				if(i == idx)
				{	Y[i] = new.X; Y[(i+1):n] = X[i:length(X)]; flag = 1; }
				else
				{	
					if (flag == 0)
					{	Y[i] = X[i]	}
				}
			}
		}
		return(Y)
	}
	else
	{	return(X)	}
}

insert.in.vector2 = function(X,new.X, idx)
{
	if  ( (length(idx) == 0) | (length(new.X) == 0) )
	{	return(X)	}
	else
	{
		Y = vector(length = (length(X) + length(idx)) )
		n = length(Y)
		for (j in 1:length(idx))
		{	Y[idx[j]] = new.X[j] 	}
		Y[-idx] = X
			
		return(Y)
	}
}


#count values 
count.factor <- function(X) 
{
	Xtable = table(X)
	values  = names(Xtable)
	row.names(Xtable) = NULL
	nb = Xtable
	return(list(values = values, numbers = nb))
}

# timer
timer = function(f)
{
	T1 = proc.time()#Sys.time()
	object = f
	T2 = proc.time()
	return(list(time.elapse = T2-T1, object = object))
}

filter.object = function(X)
{
	if (is.list(X))
	{
		if (!is.null(X$time.elapse)) {  return(X$object)	}
		else { return(X) }
	}
	else
	{ return(X)	}
}

filter.forest = function(X)
{
	if (is.list(X))
	{
		if (!is.null(X$time.elapse)) 
		{  
			if (!is.null(X$object$forest)) 	{	return(X$object$forest)	}
			else { (X$object) }
		}
		else 
		{ 
			if (!is.null(X$forest))
			{	return(X$forest) }
			else
			{	return(X)	}
		}
	}
	else
	{ return(X)	}
}

# standardize between [0,1]  or [-1,1]
standardize = function(X)
{
	if (is.vector(X))
	{	Y = standardize_vect(X,sign(min(X)  +0.00001))	}
	else
	{
		L = nrow(X); C = ncol(X);
		Y = matrix( data = NA, nrow = L , ncol = C)
		for (i in 1:C)
		{	Y[,i] =  standardize_vect(X[,i],sign(min(X[,i]) +0.00001))	}
		
		if (!is.null((colnames(X)[1])))
		{ colnames(Y) = colnames(X)	}
	}
	return(Y)
}

standardize_vect = function(X,neg)
{
	if (length(unique(X))  > 1)
	{
		n = length(X);
		max_X = max(X)
		min_X = min(X)
		Y = (X - min_X )/(max_X - min_X );
		if (neg == -1)
		{
			pos_neg = which (X < 0 )
			if ( sum(pos_neg) > 0)
			{	Y[pos_neg] = ( -X[pos_neg] + max_X )/( min_X - max_X );	}
			Y = Y[1:n];
		}
		return(Y)
	}
	else
	{ return(X) }
}


# init random train/test samples of any size
init_values <- function(X, Y = NULL, sample.size = 0.5, data.splitting = "ALL", unit.scaling = FALSE, scaling = FALSE, regression = FALSE)
{
	set.seed(sample(10000,1))
	if (is.vector(X)) {  n = length(X) } else { n = nrow(X) }
	indices = define_train_test_sets(n,sample.size,data.splitting)
	
	if (unit.scaling)
	{	X = standardize(X)	}
	
	if (scaling)
	{   X = scale(X); if (!is.null(Y) & regression) { Y = scale(Y) }	}
	
	
	if (data.splitting == "train")
	{
		xtrain = if (is.vector(X)) { X[indices$train] } else { X[indices$train,] }
		if (!is.null(Y))
		{ 	ytrain = Y[indices$train]	}
		else
		{  ytrain = Y	}
		
		return( list(xtrain = xtrain, ytrain = ytrain, train_idx = indices$train ) )
	}
	else
	{
		if (data.splitting == "test")
		{
			xtest =if (is.vector(X)) { X[indices$test] } else { X[indices$test,] }
			if (!is.null(Y))
			{ 	ytest = Y[indices$test]	}
			else
			{   ytest = Y	}
					
			return( list( xtest = xtest, ytest = ytest, test_idx = indices$test ) )
		}
		else
		{
			xtrain = if (is.vector(X)) { X[indices$train] } else { X[indices$train,] }
			xtest = if (is.vector(X)) { X[indices$test] } else { X[indices$test,] }
			
			if (!is.null(Y))
			{	ytrain = Y[indices$train]; 	ytest = Y[indices$test]	}
			else
			{   ytrain = ytest = Y	}
						
			return( list(xtrain = xtrain, ytrain = ytrain, xtest = xtest, ytest = ytest, train_idx = indices$train, 
				test_idx = indices$test) )
		}
	}
}

define_train_test_sets = function(N,sample.size,data.splitting)
{
	set.seed(sample(100,1))
	if (data.splitting == "train")
	{	
		v = sample(N)
		Ntrain = floor(N*sample.size)
		train_indices = v[1:Ntrain]
	
		return( list(train=train_indices) )
	}
	else
	{
		v = sample(N)
		Ntrain = floor(N*sample.size)
		train_indices = v[1:Ntrain]
		test_indices = v[(Ntrain+1):N]
	
		return( list(train = train_indices, test = test_indices) )
	}
}

extractYFromData <- function(XY, whichColForY = ncol(XY) )
{
	if (is.character(whichColForY))
	{	return(list( X = XY[,-which(colnames(XY) == whichColForY)], Y = XY[,whichColForY])) }
	else
	{   return(list( X = XY[,-whichColForY], Y = XY[,whichColForY])) }
}

getOddEven = function(X, div = 2)
{
	even = which( (X%%div) == 0)
	odd = which(  (X%%div) == 1)
	
	if (length(even) > 0)
	{	X.even = X[even]	}
	else
	{	X.even = 0	}
	
	if (length(odd) > 0)
	{	X.odd = X[odd]	}
	else
	{	X.odd = 0	}
	
	return(list(even = X.even, odd = X.odd))
}

min_or_max = function(X)
{
	s =sample (c(0,1),1)
	if (s)
	{ return(min(X))	}
	else
	{ return(max(X))	}

}

# Correlations 
getCorr = function(X, mid = 1/3, high= 2/3)
{
	if ( high <= mid)
	{  mid = 0.5*high	}
	
	zeros =vector()
	for (j in 1:ncol(X))
	{ zeros[j] = length(which(X[,j] == 0))/nrow(X)	}
	
	corr.matrix = cor(X)
	diag(corr.matrix) = 0
		
	high_correl = which( abs(corr.matrix) >= high, arr.ind = TRUE )
	mid_correl = which( (abs(corr.matrix) >= mid) & (abs(corr.matrix) < high), arr.ind = TRUE )
	low_correl = which( abs(corr.matrix) < mid, arr.ind = TRUE )
	
		
	if  (dim(high_correl)[1] > 0)  
	{
		high.correl.values.idx = find.first.idx(corr.matrix[high_correl])
		high.correl.values = corr.matrix[high_correl]
		var_name = cbind(high_correl, round(high.correl.values,4))[high.correl.values.idx ,]
		
		if (is.vector(var_name))
		{ var_name = as.matrix(t(var_name))	}
		
		var_name2 = var_name[match(sort(unique(var_name[,2])), var_name[,2]),]
		
		if (is.vector(var_name2))
		{ var_name2 = as.matrix(t(var_name2))	}
		
		if (!is.null(colnames(X))) 
		{	var_name = cbind(colnames(X)[var_name2[,1]], colnames(X)[var_name2[,2]], var_name2[,3])		}
		else
		{  var_name = var_name2	}
	}
	else
	{	var_name = var_name2 = paste("no correlation > ",high)	} 
		
	
	score = list(corr.matrix = corr.matrix, high.correl = high_correl, mid.correl = mid_correl, low.correl = low_correl, 
		high.correl.idx = var_name2, high.correl.var.names = var_name, zeros = round(zeros,6));
		
	return(score)
}


rm.correlation = function(X, corr.threshold = 2/3, zeros.threshold = 0.45, repeated.var = FALSE, suppress = TRUE)
{
   
    correlation.levels = getCorr(X, high = corr.threshold) 
      
    if (is.character(correlation.levels$high.correl.idx)) 
	{	return("lower correlation for all variables")	}
	else
	{	
		L1 = correlation.levels$zeros[ correlation.levels$high.correl.idx[,1]]
		L2 = correlation.levels$zeros[ correlation.levels$high.correl.idx[,2]]
		
		zeros.idx = which( ((L1+L2)/2 < zeros.threshold) | ((L1+L2) > (1+zeros.threshold)) )
		
		if ( length(zeros.idx) > 0)
		{	
			rm.idx = count.factor(correlation.levels$high.correl.idx[zeros.idx,1])
			
			if (repeated.var)
			{	rm.idx = as.numeric(rm.idx$values[which(rm.idx$numbers > 1)])	}
			else
			{ 	rm.idx = as.numeric(rm.idx$values)	}	
									
			if (is.null(colnames(X)))
			{	var_name = "no colnames"	}
			else
			{ 	var_name = colnames(X)[rm.idx]  }
					
			if (suppress)
			{	X = X[,-rm.idx]  }	
			
			high.corr.matrix = correlation.levels$high.correl.var.names	
		
			object.list = list(X = X, var_name = var_name, corr.var = rm.idx, high.corr.matrix = high.corr.matrix)
		}
		else
		{ 	object.list = "no available true correlation because of zeros. Try increase threshold"	}
		
		return(object.list)		
	}
}

find.root = function(f,a,b,...)
{
	while ((b-a) > 0.0001 )
	{
		c = (a+b)/2
		if ((f(b,...)*f(c,...)) <= 0)
		{	a = c }
		else
		{ b = c }
	}
	return(list(a=a, b=b))
}

Id <- function(X) X 

# Treatment of categorical variables
dummy.recode = function(X, which.col, threshold = 0.05, allValues = TRUE)
{
	if (is.numeric(which.col))
	{	names.X = colnames(X)[which.col]	}
	else
	{  names.X = which.col	}
	
	X = X[,which.col]
	Y = NULL
	n = length(X)
	cf.X = count.factor(X)
	cat.values = as.numeric(cf.X$values)
	cat.numbers =  as.numeric(cf.X$numbers)
	
	if (length(cat.numbers) == 1)
	{  X = as.matrix(X); colnames(X)[1] = names.X; return(X) }
	else
	{		
		merge.values.idx = which(cat.numbers <  threshold*n)
		
		if (length(merge.values.idx) > 0)
		{ 	
			cat.values[merge.values.idx]  =  cat.values[merge.values.idx[1]] 
			cat.values = unique(cat.values)
			tmp.cat.numbers =  sum(cat.numbers[merge.values.idx])
			cat.numbers[merge.values.idx] = 1
			cat.numbers = unique(cat.numbers)
			cat.numbers[which(cat.numbers == 1)] = tmp.cat.numbers
		}
				
		dummy.l = length(cat.values) + (allValues - 1)
		tmp.X = vector(length = length(X))
		tmp.names = names.X 
		
		if (dummy.l == 0)
		{	X = rep(1, length(X));  X = as.matrix(X); colnames(X)[1] = names.X; return(X)  }
		else
		{
			for (i in 1:(dummy.l))	
			{
				idx = which(X == cat.values[i] )
				if (length(merge.values.idx) > 0)
				{
					if (cat.values[i] == merge.values.idx[1])
					{ 
						for (j in 2:length(merge.values.idx))
						{	idx = c(idx, which(X == merge.values.idx[j]))	}
					}		
				}
				tmp.X[idx] =  1
				tmp.X[-idx] = 0		
				Y = cbind(Y,tmp.X)	
				tmp.names[i] = paste(names.X,".", cat.values[i],sep = "")
			}	

			colnames(Y) = tmp.names
			return(Y)
		}		
	}
}

inDummies <- function(X, dummies, inDummiesThreshold = 0.05, inDummiesAllValues = TRUE )
{
	X.dummies = NULL
	for (j in 1:length(dummies)) {	X.dummies = cbind(X.dummies, dummy.recode(X, dummies[j], threshold = inDummiesThreshold, allValues = inDummiesAllValues)) }
	
	if (is.numeric(dummies[1])) {  return(cbind(X[,-dummies], X.dummies)) }
	else {   return(cbind(X[,-match(dummies, colnames(X))], X.dummies))	}
}

rmNoise =  function(X, catVariables = NULL, threshold = 0.5, inGame = FALSE, columnOnly = TRUE)
{
	n = nrow(X)
	if (is.null(catVariables)) {   catVariables = 1:ncol(X)	}
	
	if (inGame)
	{
		for (j in catVariables)
		{
			idxInf = which( table(X[,catVariables[j]])/n < threshold )
			if (length(idxInf) > 0)
			{ 
				for (i in 1:length(idxInf))
				{	
					rmIdx = which(X[,j] == as.numeric(names(idxInf[i])))
					if (length(rmIdx) > 0)	{ 	X = X[-rmIdx,]	}
				}
			}			
			n = nrow(X)
		}
	}
	else
	{
		rmIdxAll = NULL
		for (j in catVariables)
		{
			idxInf = which( table(X[,catVariables[j]])/n > threshold )
			if (length(idxInf) > 0)
			{ 
				if (columnOnly)
				{	rmIdxAll = c(rmIdxAll, catVariables[j])	  }
				else
				{
				
					for (i in 1:length(idxInf))
					{	
						rmIdx = which(X[,j] == as.numeric(names(idxInf[i])))
						if (length(rmIdx) > 0)	{ 	rmIdxAll = c(rmIdxAll, rmIdx)	}
					}
				}
			}			
		}
		
		if (!is.null(rmIdxAll))	{ rmIdxAll = unique(rmIdxAll)	}
		if (columnOnly)	{  X = X[,-rmIdxAll] }
		else {	X = X[-rmIdxAll,]	}
	}
	
	return(X)
}


# category2Quantile <- function(X, nQuantile = 20)
# {
	# Z = X
	
	# for (j in 1:ncol(X))
	# {
		# limQ = min(X[,j])
		# for (i in 1:nQuantile)
		# {	
			# Q = quantile(X[,j], i/nQuantile)
			# Z[which(Z[,j] <= Q & Z[,j] >= limQ), j] = i-1
			# limQ = Q	
		# }
	
	# }
	# return(Z)
# }

# category2Proba <- function(X,Y, whichColumns = NULL, regression = FALSE, threads = "auto")
# {
	# if (is.null(whichColumns)) {  whichColumns = 1:ncol(X) }
	# p = length(whichColumns)
	# n = nrow(X)
	# ghostIdx = 1:n
	# Z = X
	
	# for (j in  whichColumns)
	# {
		# XX = cbind(X[,j],ghostIdx)
		# NAIdx = which.is.na(X[,j])
		# if (length(NAIdx) > 0)
		# {	
			# catValues = rmNA(unique(X[-NAIdx,j])) 				
			# for (i in seq_along(catValues))
			# {
				# idx = which(X[-NAIdx,j] == catValues[i])
				# trueIdx = XX[-NAIdx,j][idx]
				# if (!regression)
				# {	Z[trueIdx, j] = max(table(Y[trueIdx]))/length(idx)	}
				# else
				# {	Z[trueIdx, j] = sum(Y[trueIdx])/length(idx) }
			# }
		# }
		# else
		# {
			# catValues = unique(X[,j])
			# for (i in seq_along(catValues))
			# {
				# idx = which(X[,j] == catValues[i])
				# if (!regression)
				# { 	Z[idx, j]= max(table(Y[idx]))/length(idx)	} #Z[idx ,j]
				# else
				# {	Z[idx,j] = sum(Y[idx])/length(idx) }
			# }
		# }
	# }
	# return(Z)
# }
		
# imputeCategoryForTestData <- function(Xtrain.imputed, Xtrain.original, Xtest, impute = TRUE)
# {
	# Xtest.imputed = matrix(NA, ncol = ncol(Xtrain.imputed), nrow = nrow(Xtest))
	# for (j in 1:ncol(Xtrain.original))
	# {
		# uniqueCat = unique(Xtrain.original[,j])
		# for( i in 1:length(uniqueCat) )
		# {
			# uniqueIdx = which(Xtest[,j] == uniqueCat[i])
			
			# if (length(uniqueIdx) > 0)
			# {
				# uniqueIdx2 = which(Xtrain.original[,j] == uniqueCat[i])[1]
				# Xtest.imputed[uniqueIdx,j] = Xtrain.imputed[uniqueIdx2,j]
			# }
		# }
	# }
	
	# if (impute)
	# {	
		# if (length(which.is.na(Xtest.imputed)) > 0)
		# {	Xtest.imputed = na.impute(Xtest.imputed) }
	# }
		
	# return(Xtest.imputed)
# }
			
# categoryCombination <- function(X, Y, idxBegin = 1, idxtoRemove = NULL)
# {
	# seqCol = (1:ncol(X))[-c(idxBegin, idxtoRemove)]
	# Z = matrix(data = 0, ncol = length(seqCol), nrow = nrow(X))
	# l = 1
	# for (j in seqCol)
	# {
		# matchSameCat = which(X[,idxBegin] == X[,j])
		# if (length(matchSameCat) > 0)
		# {
			# Z[matchSameCat, l] = max(table(rmNA(Y[matchSameCat])))/length(matchSameCat)
			# l = l + 1
		# }
	# }
	
	# return(Z)
# }
		
permuteCatValues = function(X, catVariables)
{
	if (is.vector(X))
	{ X = as.matrix(X); catVariables = 1 }

	XX = X 
		
	for (j  in catVariables)
	{
		tmpX = Z = X[,j]
		XTable  = count.factor(Z)
		XCat = as.numeric(XTable$values)
		XNb = XTable$numbers
		
		nRepCat = length(XCat)
		randomCat = sample(XCat, nRepCat)
		
		for (i in 1:nRepCat)
		{
			if (XNb[i] > 0)
			{
				oldIdx = which(Z == XCat[i])
				tmpX[oldIdx] = randomCat[i]
			}
		}
		XX[,j] = tmpX
	}
	
	return(XX)
}

generic.log <- function(X, all.log = FALSE)
{
	if (all.log)
	{
		for (j in 1:ncol(X)) { 	X[,j] = log(X[,j] + abs(min(X[,j])) +1) }
	}
	else
	{
		for (j in 1:ncol(X)) { 	X[,j] = if (min(X[,j]) >= 0) { log(X[,j] + 1) } else { X[,j] } }
	}
	return(X)
}

smoothing.log = function(X, filtering = median, k = 1)
{
	zeros = which(X == 0)
	if (length(zeros) > 0)
	{
		min.X = filtering(X[-zeros])
		replace.zeros = runif(length(zeros),0, k*min.X)
		X[zeros] = replace.zeros 
	}
	
	return(X)
	
}

generic.smoothing.log = function(X, threshold = 0.05, filtering = median ,  k = 1)
{
	if (is.vector(X) )
	{	n  = length(X);  p = 1; X =	as.matrix(X)	}
	else
	{ 	p = ncol(X); n = nrow(X) }
	
	Z = X
	
	for (j in 1:p)
	{ 
		min_X = min(X[,j])
		if ( (max(X[,j]) > 0) & ( min_X < 0) )
		{	
			neg.idx = which(X[,j] < 0)
			pos.idx = which(X[,j] < 0)
			
			if (length(neg.idx) <= threshold*n)
			{	X[neg.idx,j] =	0;  Z[,j] = log(smoothing.log(abs(X[,j]), filtering = median ,  k = 1)) 	}
		
			if (length(pos.idx) <= threshold*n)
			{	X[pos.idx,j] =	0;  Z[,j] = log(smoothing.log(abs(X[,j]), filtering = median ,  k = 1)) 	}
		}
		else
		{	
			
			if ( (min(X[,j]) != 0) & (max(X[,j]) != 0) )
			{			
				if ( (min(X[,j]) == 0) | (max(X[,j]) == 0) )
				{	Z[,j] = log(smoothing.log(abs(X[,j]), filtering = median,  k = 1))	}
			}
		}
	}
	
	return(Z)
}


confusion.matrix = function(predY, Y)
{
	if ( (is.vector(predY) & !is.vector(Y)) | (!is.vector(predY) & is.vector(Y)) ) 
	stop("Predictions and responses do not have the same class.")
	if (length(predY) != length(Y)) stop("Predictions and responses do not have the same length.")
	if ( (class(Y) == "factor") | (class(predY) == "factor") ) 
	{  
		trueClasses = levels(as.factor(Y))
		predY = as.numeric(predY);  Y = as.numeric(Y)
		classes = 1:length(trueClasses)
	}
	else
	{
		trueClasses = 1:length(tabulate(Y))
		predClasses = 1:length(tabulate(predY))
		if (length(predClasses) > length(trueClasses)) 	{  trueClasses = predClasses }
		classes = trueClasses
	}
	
	n = length(Y)
	k = length(classes)
	TP = vector(length = length(classes))
		
	if (k == 2)
	{ 	FP = vector(length = length(classes))	}
	else
	{ 	FP = matrix( data = 0, ncol = k, nrow = k ) } 
	
	for (i in 1:k)
	{
		TP[i] = length(which((predY == classes[i]) & (Y == classes[i])))
	   
	    if (k > 2)
	    {
			out.classes = which(classes != i)
			FP[i,i] = TP[i]
			for (j in 1:(k-1))
			{	  FP[i,out.classes[j]] =  length(which((predY == classes[i]) & (Y == classes[out.classes[j]])))	}
		}
		else
		{	FP[i] = length(which((predY == classes[i]) & (Y !=classes[i])))		}
		
	}
		
	if (k == 2)
	{	ConfMat = matrix(data = cbind(TP, FP) , nrow = k, ncol = k); ConfMat[2,] = t(c(FP[2],TP[2])) 		}
	else
	{	ConfMat = FP	}
	
	rownames(ConfMat) = trueClasses
	colnames(ConfMat) = trueClasses
	class.error = round(1- diag(ConfMat)/rowSums(ConfMat),4)
	
	ConfMat = cbind(ConfMat,class.error)
	names(dimnames(ConfMat)) = c("Prediction", "Reference")

	return(ConfMat)
}


generalization.error = function(conf.matrix) 1 - sum(diag(conf.matrix[,-ncol((conf.matrix))]))/sum(conf.matrix[,-ncol((conf.matrix))])

overSampling = function(Y, which.class = 1, proportion = 0.5)
{
	n = length(Y)
	S1 <- find.idx(which.class, Y, sorting = FALSE)
	if (proportion == -1)
	{	S1.newSize = length(Y[-S1])	}
	else
	{	S1.newSize = floor((1 + proportion)*length(S1)) }
	
	S1.sample = sample(S1, S1.newSize, replace = ifelse(length(S1) > S1.newSize, FALSE, TRUE) ) 
	S2 = (1:n)[-S1]
	if (proportion == -1)
	{	S2.sample = S2 }
	else
	{ 	S2.sample = sample(S2, length(S2) + length(S1) - S1.newSize, replace = TRUE)  }#
			
	return( list(C1 = sort(S1.sample), C2 = sort(S2.sample)) )
}

outputPerturbationSampling = function(Y, whichClass = 1, sampleSize = 1/4, regression = FALSE)
{
	n = length(Y)
	
	if (regression)
	{
		randomIdx =  sample(1:n, min(n, floor(sampleSize*n) + 1), replace = TRUE)
		meanY <- sum(Y[randomIdx])/length(Y[randomIdx])
		minY <- min(Y)
		maxY <- max(Y)
		
		sdY <- if (sum(abs(Y[randomIdx])) == 0) { 0.01 } else { max(0.01, sd(Y[randomIdx])) }
			
			
		randomData = rnorm(length(randomIdx), meanY, 3*sdY)
		if ( (minY < 0) & (maxY > 0)) 
		{	Y[randomIdx] =  randomData  }
				
		if ( (minY >= 0) & (maxY > 0))
		{	Y[randomIdx] =  abs(randomData) } 
		
		if ( (minY < 0 ) & (maxY <= 0))
		{	Y[randomIdx] =  -abs(randomData) } 		
	}
	else
	{	
		if (whichClass == -1)	{ 	whichClass = 1	}
		
		idxClass = which(Y == whichClass)
		
		perturbClassidx = sample(idxClass, floor(min(sampleSize,0.5)*length(idxClass)))
		perturbClassY  = sample(Y[which(Y !=  whichClass)], length(perturbClassidx), replace = TRUE)
		Y[perturbClassidx] = perturbClassY
	}
	
	return(Y)
}


keep.index <- function(X,idx) if (!is.null(dim(X))) { (1:nrow(X))[idx] } else  { (1:length(X))[idx] }

roc.curve <- function(X, Y, classes, positive = classes[2], ranking.threshold = 0, ranking.values = 0, falseDiscoveryRate = FALSE, plotting = TRUE,  printValues = TRUE, Beta = 1)
{
	#require(pROC)
	par(bg = "white")
	classes = sort(classes)
	
	if (length(classes) != 2)
	{ stop("Please provide two classes for plotting ROC curve.") }
	
	if (is.vector(X) & (!is.numeric(X)) & !(is.list(X)) )
	{	X <- vector2factor(X) }
	
	if (is.factor(X))
	{ 	X <- factor2vector(X)$vector	}
	
	if (is.list(X))
	{
		if (class(X) == "randomUniformForest")
		{	
			if (!is.null(X$predictionObject))
			{	X <- X$predictionObject	}
			else
			{   X <- X$forest }
		}		
	}
		
	YObject = factor2vector(vector2factor(Y))
	
	Y = YObject$vector
	newClasses = rep(0,2)
	newClasses[1] = find.idx(classes[1], YObject$factors)
	newClasses[2] = find.idx(classes[2], YObject$factors)
	positive = newClasses[2]
	classes = newClasses	
	
	if (falseDiscoveryRate)
	{		
		if (is.list(X))
		{
			class.object = if (!is.null(X$all.votes)) { majorityClass(X$all.votes, classes) } 
			else { majorityClass(X$OOB.votes, classes) }  
			pred.classes = if (!is.null(X$majority.vote)) { X$majority.vote } else { X$OOB.predicts }
			class.probas = (class.object$class.counts/rowSums(class.object$class.counts))[,1]
		}
		else
		{  
			class.probas = 0.5; 
			pred.classes = X 
			
			if (length(X) != length(Y))
			{ stop("X is not a vector, or length of X and Y are not the same.\n") }
		}
	}
	else
	{ 	
		pred.classes = if (is.numeric(X)) { X } 
		else 
		{  
			if (!is.null(X$all.votes))	{ X$majority.vote }
			else  {	X$OOB.predicts }	
		}
		class.probas = ranking.values;   
		ranking.threshold = 0 
	}
		
	pos.idx = which(pred.classes == positive & class.probas >=  ranking.threshold )
	n_pos.idx = length(pos.idx)
	
	if (n_pos.idx == 0)
	{
		cat("\nNo positive case found.\n")
		return(list(cost = 0, threshold = ranking.threshold, accuracy = 0, expected.matches = 0))
	}
	else
	{
		conf.mat = confusion.matrix(pred.classes,Y)
		n2 = diag(conf.mat[,-ncol(conf.mat)])[nrow(conf.mat)] # true positives
		n1 = length(pos.idx) - n2		# false positives
		n3 = conf.mat[1,2]	# false negatives
		n4 = conf.mat[1,1]   # true negatives
		
		inv_n23 = 1/(n2 + n3)
		inv_n12 = 1/(n1 + n2)
		if (falseDiscoveryRate)
		{
			VP = array(0, n_pos.idx+1)
			FP = array(1, n_pos.idx+1)

			for(i in 1:n_pos.idx)
			{
				if (pred.classes[pos.idx[i]] == Y[pos.idx[i]])	{ VP[i+1] = VP[i] + inv_n23; FP[i+1] = FP[i] }
				else  {	FP[i+1] = FP[i] - inv_n12;  VP[i+1] = VP[i] }
			}
			FP.new = c(FP,0)
			VP.new = c(VP,1)
			
			n_pos = length(VP.new)
			aupr = rep(0,n_pos)
			for (j in 2:n_pos)
			{	
				DV = (VP.new[j] - VP.new[j-1])
				DF = (FP.new[j-1] - FP.new[j] )
				#DF = (1-FP.new[j] - (1-FP.new[j-1]))
				aupr[j] = { (1 -  VP.new[j])*DF + 0.5*DV*DF }
			}
			aupr = 1- sum(aupr)
		}
		else
		{
			# for ROC curve : VP rate = Sensitivity = VP/(VP + FN) ;  specificity = VN/(VN + FP) = 1 - FP rate 
			# and Sensitivity = True positives among all; Specificity = True negatives among all
			VP = FP = array(0, n_pos.idx+1)
			inv_n14 = 1/(n1 + n4)
			for(i in 1:n_pos.idx)
			{
				if(pred.classes[pos.idx[i]] == Y[pos.idx[i]] )
				{	VP[i+1] = VP[i] + inv_n23; FP[i+1] = FP[i]	}
				else
				{	FP[i+1] = FP[i] + inv_n14;  VP[i+1] = VP[i] }
			}
			FP.new = c(FP,1)
			VP.new = c(VP,1)
			AUC.est = round(pROC::auc(Y, pred.classes ),4)
		}
			
		if (plotting)
		{
			if (ranking.threshold == 0)
			{	
				F1.est = fScore(conf.mat, Beta = Beta)
				
				if(!falseDiscoveryRate)
				{
					plot(FP.new, VP.new, xlab = "False Positive rate [1 - Specificity]", ylab = "Sensitivity [True Positive rate]", 
					type='l', col = sample(1:10,1), lwd = 2)
					title(main = "ROC curve", sub = paste("AUC: ",  AUC.est, ". Sensitivity: ", round(100*VP.new[length(pos.idx)+1],2), "%",
					". False Positive rate: ", round(100*FP.new[length(pos.idx)+1],2), "%", sep=""), col.sub = "blue", cex.sub = 0.85)
					points(c(0,1), c(0,1), type='l', col='grey')
					abline(v = FP.new[length(pos.idx)+1], col='purple', lty = 3)
					abline(h = VP.new[length(pos.idx)+1],col='purple', lty = 3)
					if (printValues) { cat("AUC: ",  AUC.est, "\n") }
				}
				else
				{
					plot(VP.new, FP.new, xlab = "Sensitivity [True Positive rate]", ylab = "Precision [1 - False Discovery rate ]", 
					type='l', col = sample(1:10,1), lwd = 2)
					title(main = "Precision-Recall curve", sub = paste("\nProbability of positive class", " > ", 0.5, ". Precision: ", 
					round(100*FP.new[length(pos.idx)+1],2), "%", ". Sensitivity: ", round(100*VP.new[length(pos.idx)+1],2), "%",  ". F", Beta, "-Score: ", round(F1.est,4), ". AUPR: ", round(aupr, 4), sep = ""), col.sub = "blue", cex.sub = 0.85)
					points(c(0,1), c(1,0), type='l', col='grey')
					abline(v = VP.new[length(pos.idx)+1], col='purple', lty = 3)
					abline(h = FP.new[length(pos.idx)+1], col='purple', lty = 3)
					if (printValues) 
					{
						cat("AUPR (Area Under the Precision-Recall curve) :", round(aupr,4), "\n")
						cat("F1 score: ", round(F1.est,4) , "\n")
					}
				}
				#grid()
				if (printValues) { print(conf.mat) }
			}
			else
			{	
				pred.classes[-pos.idx] = classes[-which(classes == positive)]
				conf.mat = confusion.matrix(pred.classes,Y)
				F1.est = fScore(conf.mat)
				if (printValues) 
				{
					cat("F1 ", F1.est, "\n")
					print(conf.mat)
				}
				plot(VP.new, FP.new, xlab = "Sensitivity [1 - True positives Rate]", ylab = "Precision [1 - False Discovery Rate]", main = paste("Precision-recall curve ", "[Probability of positive class", " > ", ranking.threshold, "]", sep=""), type='l', col = sample(1:10,1), lwd = 2)
				points(c(0,1), c(1,0), type='l', col='grey')
				abline(v = VP.new[length(pos.idx)+1], col='purple', lty = 3)
				abline(h = FP.new[length(pos.idx)+1], col='purple', lty = 3)
				#grid()
			}
		}
		else
		{
			return(list(cost = round((1-(n1+n2)/sum(pred.classes == positive))*100,2), threshold =  ranking.threshold, 
				accuracy = round(( 1 - (1 - F1.est)/conf.mat[nrow(conf.mat), ncol(conf.mat)])*100,2), 
				expected.matches = round(((n1 + n2)/sum(Y == positive))*100,2)))
		}
	}
}


myAUC <- function(pred.classes, Y, positive = 2, falseDiscoveryRate = FALSE)
{	
	pos.idx = which(pred.classes == positive)
	n_pos.idx = length(pos.idx)
	if(n_pos.idx == 0)
	{	
		cat("\nNo positive case found.\n")
		return(list(auc = 0, VP = 0, FP = 0))
	}
	else
	{
		conf.mat = confusion.matrix(pred.classes,Y)
		
		n2 = diag(conf.mat[,-ncol(conf.mat)])[nrow(conf.mat)] # true positives
		n1 = n_pos.idx - n2	# false positives
		n3 = conf.mat[1,2]	# false negatives
		n4 = conf.mat[1,1]  # true negatives
		
		inv_n23 = 1/(n2 + n3)
		inv_n12 = 1/(n1 + n2)
		#compute Area under Precision-Recall Curve
		if (falseDiscoveryRate) 
		{
			VP = array(0, n_pos.idx+1)
			FP = array(1, n_pos.idx+1)
			for(i in 1:n_pos.idx)
			{
				if (pred.classes[pos.idx[i]] == Y[pos.idx[i]] )  { VP[i+1] = VP[i] + inv_n23; FP[i+1] = FP[i]	}
				else {	FP[i+1] = FP[i] - inv_n12;  VP[i+1] = VP[i] }
			}
			FP.new = c(FP,0)
			VP.new = c(VP,1)
			
			n_pos = length(VP.new)
			auc = rep(0,n_pos)
			for (j in 2:n_pos)
			{	
				DV = (VP.new[j] - VP.new[j-1])
				DF = (FP.new[j-1] - FP.new[j] )
				#DF = (1-FP.new[j] - (1-FP.new[j-1]))
				auc[j] = { (1 -  VP.new[j])*DF + 0.5*DV*DF }
			}
			auc = 1- sum(auc)
		}
		# compute Area under ROC Curve
		else
		{
			# for ROC curve : VP rate = Sensitivity = VP/(VP + FN);  specificity = VN/(VN + FP) = 1 - FP rate 
			# and Sensitivity = True positives among all; Specificity = True negatives among all
			# rOC plots True positive rate (y-axis) vs false positive rate (x-axis)
			VP = FP = array(0, length(pos.idx)+1)
			inv_n14 = 1/(n1 + n4)
			for (i in 1:length(pos.idx))
			{
				if (pred.classes[pos.idx[i]] == Y[pos.idx[i]] ) {	VP[i+1] = VP[i] + inv_n23; FP[i+1] = FP[i]	}
				else {	FP[i+1] = FP[i] + inv_n14; VP[i+1] = VP[i] }
			}
			FP.new = c(FP,1)
			VP.new = c(VP,1)
			
			n_pos = length(VP.new)
			auc = rep(0,n_pos)
			for (j in 2:n_pos)
			{	
				DV = (VP.new[j] - VP.new[j-1])
				DF = (FP.new[j] - FP.new[j-1])
				auc[j] = (1 - VP.new[j])*DF + 0.5*DV*DF
			}
			auc = 1- sum(auc)
		}
		
		return(list(auc = auc, VP = VP.new, FP = FP.new))
	}
}


optimizeFalsePositives = function(X, Y, classes, o.positive = classes[2], o.ranking.values = 0, o.falseDiscoveryRate = TRUE, stepping = 0.01)
{
	if (o.falseDiscoveryRate)
	{ 	threshold = seq(0.51, 0.9, by = stepping)	}
	else
	{	 threshold = seq( median(o.ranking.values), max(o.ranking.values), by = 0.01)	}
	
	tmp.AUC = 0.01
	cost = accuracy = expected.matches = AUC = vector()
	
	i = 1; n = length(threshold)
	while( (tmp.AUC < 1) & (tmp.AUC > 0) & (i <= n))
	{
		ROC.object = roc.curve(X, Y, classes, positive = o.positive, ranking.threshold = threshold[i], ranking.values = o.ranking.values, falseDiscoveryRate = o.falseDiscoveryRate, plotting = FALSE)
		
		cost[i] = ROC.object$cost 
		accuracy[i] = ROC.object$accuracy 
		expected.matches[i] = ROC.object$expected.matches
		AUC[i] = ROC.object$AUC
		tmp.AUC = AUC[i]
		i = i + 1
	}
	threshold = threshold[1:(i-1)]
	
	plot(threshold, cost, type='l')
	points(threshold, accuracy, col='red', type='l')
	points(threshold, AUC*100, col='green', type='l')
	points(threshold, expected.matches, col='blue', type='l')
	grid()

	return(list(cost = cost, threshold = threshold,  accuracy = accuracy , expected.matches = expected.matches,  AUC = AUC))
}

someErrorType <- function(modelPrediction, Y, regression = FALSE)
{
	#require(pROC)	
	if (regression)
	{  
		n = length(Y)
		MSE <- L2Dist(Y, modelPrediction)/n
		return(list(error = MSE, meanAbsResiduals = sum(abs(modelPrediction - Y))/n, Residuals = summary(modelPrediction - Y), percent.varExplained = max(0, 100*round(1 - MSE/var(Y),4)) ) )
	}
	else
	{  
		confusion = confusion.matrix(modelPrediction,Y)
		if (length(unique(Y)) == 2 )
		{	
			return(list( error = generalization.error(confusion), confusion = confusion, 
			AUC = pROC::auc(Y, as.numeric(modelPrediction)), AUPR = myAUC(as.numeric(modelPrediction), as.numeric(Y), falseDiscoveryRate = TRUE)$auc) )
		}
		else
		{	return(list( error = generalization.error(confusion), confusion = confusion))	}
	}	
}

gMean <- function(confusionMatrix, precision = FALSE)  
{
	p = ncol(confusionMatrix)
	if (precision) 	{  acc = rmNA(diag(confusionMatrix[,-p])/rowSums(confusionMatrix[,-p])) }	
	else { 	acc = rmNA(diag(confusionMatrix[,-p])/colSums(confusionMatrix[,-p]))  }
	idx = which(acc == 0)
	if (length(idx) > 0) { acc = acc[-idx]	}
	nbAcc = length(acc)
	return( (prod(acc))^(1/nbAcc))
}

fScore <- function(confusionMatrix, Beta = 1) 
{
	confusionMatrix = confusionMatrix[,-ncol(confusionMatrix)]
	precision = confusionMatrix[2,2]/(confusionMatrix[2,1] +  confusionMatrix[2,2])
	recall = confusionMatrix[2,2]/(confusionMatrix[1,2] +  confusionMatrix[2,2])
	
	return( (1 + Beta^2)*(precision*recall)/(Beta^2*precision + recall) )
}

expectedSquaredBias <- function(Y, Ypred)  (mean(Y - mean(Ypred)))^2

randomWhichMax <- function(X)
{
    set.seed(sample(1e9,1))
	Y = seq_along(X)[X == max(X)]
    if (length(Y) == 1) { Y } else  { sample(Y,1) } 
}
	
# import
estimaterequiredSampleSize <- function(accuracy, dOfAcc) -log( dOfAcc/2)/(2*(1- accuracy)^2)  
estimatePredictionAccuracy <- function(n,  dOfAcc = 0.01)
{
	virtualRoot <- function(accuracy, n = n, dOfAcc = dOfAcc)  n - estimaterequiredSampleSize(accuracy, dOfAcc)
	accuracy <- mean(unlist(find.root(virtualRoot, 0.5, 1, n = n, dOfAcc = dOfAcc)))
	
	return(accuracy)
}

subEstimaterequiredSampleSize <- function(accuracy, n) pbinom( floor(n*0.5)-floor(n*(1- accuracy)), size = floor(n), prob=0.5) 
binomialrequiredSampleSize <- function(accuracy, dOfAcc) 
{ 
  #require(gtools)
	
  r <- c(1,2*estimaterequiredSampleSize(accuracy, dOfAcc))
  v <- binsearch( function(n) { subEstimaterequiredSampleSize(accuracy, n) - dOfAcc }, range = r, lower = min(r), upper = max(r))

  return(v$where[[length(v$where)]])
}
# End of import

setManyDatasets <- function(X, n.times, dimension = FALSE, replace = FALSE, subsample = FALSE, randomCut = FALSE)
{
	X.list = vector("list", n.times)
	if (dimension > 0)
	{  
		p = ncol(X)
		minDim = sqrt(p)
		if (dimension != 1) {  toss = dimension }
		for (i in 1:(n.times))
		{	
			if (dimension == 1) { toss = sample(minDim:p,1)/p }
			X.list[[i]] = sort(init_values(1:p, sample.size = toss, data.splitting = "train")$xtrain) 
		}
	}
	else
	{
		n = nrow(X)
		if (randomCut)	{ 	idx = sample(1:n, n) }
		else { idx = 1:n }
		cuts = cut(idx, n.times, labels = FALSE)
		n.times = 1:n.times			
		if (replace & !subsample) 
		{ 
			idx = sort(sample(idx, n, replace = replace))
			for (i in seq_along(n.times))
			{	X.list[[i]] = idx[which(cuts == i)] }
		}
		else
		{
			if (subsample != 0) 
			{ 
				idx = sort(sample(idx, floor(subsample*n), replace = replace))
				for (i in seq_along(n.times))
				{	X.list[[i]] = rmNA(idx[which(cuts == i)]) }
			}
			else
			{
				for (i in seq_along(n.times))
				{	X.list[[i]] = which(cuts == i) }
			}
		}		
	}
	
	return(X.list)
}

# simulation of random gaussian matrix which params come from an uniform distribution between [-10,10]
simulationData <- function(n, p, distrib = rnorm, colinearity = FALSE)
{
	X = matrix(NA, n, p)
	for (j in 1:p)
	{
		params = runif(1, -10, 10)
		X[,j] = distrib(n, params, sqrt(abs(params)))
	}
	X <- fillVariablesNames(X)
	
	return(X)
}

difflog <- function(X) 	diff(log(X))

rollApplyFunction <- function(Y, width, by = NULL, FUN = NULL)   
{  	
    if (is.null(by)) { by = width }	
		
	n = length(Y) 
	T = min(trunc(n/width),n)
	Z = NULL
	j = 1
	
	for(i in 1:T) 
	{	
		ZTmp = Y[j:(j + width - 1)]
		idxZ = seq(1, width, by = by)
		Z = c(Z,FUN(ZTmp[idxZ]))
		j = j + width 
	} 
	
	return( na.omit(Z))
}

lagFunction <- function(X, lag = 1, FUN = NULL, inRange= FALSE)
{
	if (lag == 1) { return(FUN(X))	}
	else
	{
		idx = 1:length(X)
		lagIdx = idx + lag
       
		limitIdx = which(lagIdx == max(idx))
	   
		XIdx = cbind(lagIdx[1:limitIdx],idx[1:limitIdx])
	   
		if (inRange) {  return(apply(cbind(XIdx[,2], XIdx[,1]),1, function(Z)  FUN(X[Z[1]:Z[2]]) )) }		
	   else { return(apply(cbind(X[XIdx[,2]], X[XIdx[,1]]),1, function(Z)  FUN(Z) ))     } 
	   
	}
}

L2Dist <- function(Y,Y_)  sum( (Y - Y_)*(Y - Y_) )

L1Dist <- function(Y,Y_) sum(abs(Y - Y_))

L2.logDist <- function(Y,Y_) sum( (log(Y) - log(Y_))^2 )

HuberDist <- function(Y, Y_, a = abs(median(Y - Y_)))  sum(ifelse( abs(Y - Y_) <= a, (Y - Y_)^2, 2*a*abs(Y - Y_)- a^2))

pseudoHuberDist <- function(Y, Y_, a = abs(median(Y - Y_))) sum(a^2*(sqrt(1 + ((Y - Y_)/a)^2) - 1))

# ndcg <- function(estimatedRelevanceScoresNames, trueRelevanceScoresNames, predictionsMatrix, idx = 1:nrow(predictionsMatrix)) 
# {
	# scores = sortMatrix(predictionsMatrix[idx,], trueRelevanceScoresNames, decrease = TRUE)
	# estimatedRelevanceRanks = sortMatrix(predictionsMatrix[idx,], "rank")
	# trueRelevanceScores = scores[,trueRelevanceScoresNames]
	# scoreIdx = match(scores[,estimatedRelevanceScoresNames], estimatedRelevanceRanks[,estimatedRelevanceScoresNames]) 
	# estimatedRelevanceScores = estimatedRelevanceRanks[ scoreIdx,"majority vote"]
	# zeros = which(trueRelevanceScores == 0)
	# if (length(zeros) > 0)
	# { estimatedRelevanceScores[zeros] = 0 }
	# estimatedRelevanceScores = sortMatrix(cbind(scoreIdx, estimatedRelevanceScores),1)[,2]
  
	# limitOverEstimation = which(trueRelevanceScores > 0)
	# if (length(limitOverEstimation) > 0)
	# {	
		# estimatedRelevanceScores[scoreIdx[limitOverEstimation]] = min(estimatedRelevanceScores[scoreIdx[limitOverEstimation]], 
		# trueRelevanceScores[limitOverEstimation])
	# }
     
	# DCG <- function(y) sum( (2^y - 1)/log(1:length(y) + 1, base = 2) ) #y[1] + sum(y[-1]/log(2:length(y), base = 2))
 
	# DCGTrueScore = DCG(trueRelevanceScores)
	# DCGScore = DCG(estimatedRelevanceScores)
	# if (DCGTrueScore == 0)
	# { return (-1)  }
	# else
	# {	 return(DCGScore/DCGTrueScore) }
# }

# fullNDCG <- function(estimatedRelevanceScoresNames, trueRelevanceScoresNames, predictionsMatrix, idx = 1:nrow(predictionsMatrix), filterNoOrder = TRUE)
# {
	# if (is.null(idx))
	# { stop("Please use ndcg() function. fullNDCG needs index of ID to rank scores.") }
	
	# cases = sort(unique(idx))
	# scores = vector(length = length(cases))
	# for (i in 1:length(cases))
	# {
		# subCases = which(idx == cases[i])
		# scores[i] <- ndcg(estimatedRelevanceScoresNames, trueRelevanceScoresNames, predictionsMatrix, idx = subCases)
	# }

	# if (filterNoOrder)
	# { scores[scores == -1] == 1	}
	
	# return(scores)
# }
	
randomize = function(X)  sample(1:nrow(X), nrow(X))

# Bayesian bootstrap
# sampleDirichlet <- function(X)
# {
	# gamma.sample = rgamma(ncol(X),1)
	# dirichlet.sample = gamma.sample/sum(gamma.sample) 
	# R.x =  apply(X, 2,  function(Z) c(max(Z),min(Z)))
	
	# return(-diff(R.x)*dirichlet.sample + R.x[2,])
# }

#some colors
#require(grDevices)
# palet <- colorRampPalette(c("yellow", "red"))

bCICore <- function(X, f = mean, R = 100, conf = 0.95, largeValues = TRUE, skew = TRUE)
{
	XX = X
	bootstrapf = rep(0,R)
	X = unique(X)
	nn = length(X)
	
	replaceWithSeed <- function(f, X, n)
	{
		set.seed(sample(2*R*length(XX),1))
		f(sample(X, n, replace = FALSE))
	}
	
	bootstrapf <- replicate(R, replaceWithSeed(f, XX, nn)) 

	alpha1 = (1 - conf)/2
	alpha2 = conf + alpha1
	if (largeValues) 
	{ 
		if (!skew) {	return(c(mean(bootstrapf), quantile(bootstrapf, alpha1) + qnorm(alpha1)*sd(X)/sqrt(nn), quantile(bootstrapf, alpha2) + qnorm(alpha2)*sd(X)/sqrt(nn) )) }
		else 
		{ 
			return(c(mean(bootstrapf), quantile(bootstrapf, alpha1) - sd(X), quantile(bootstrapf, alpha2)+ sd(X)))	
		}		
	}
	else  { return(c(mean(bootstrapf), quantile(bootstrapf, alpha1), quantile(bootstrapf, alpha2) )) }
}


bCI <- function(X, f = mean, R = 100, conf = 0.95, largeValues = TRUE, skew = TRUE, threads = "auto")
{
	if (is.matrix(X)) { n = nrow(X) } else { n =length(X) }
	
	{
		#require(parallel)	
		max_threads = detectCores()
		
		if (threads == "auto")
		{	
			if (max_threads == 2) { threads = max_threads }
			else {	threads  = max(1, max_threads - 1)  }
		}
		else
		{
			if (max_threads < threads) 
			{	cat("Warning : number of threads is higher than logical threads in this computer.\n") }
		}
		
		{
			#require(doParallel)
			
			Cl = makePSOCKcluster(threads, type = "SOCK")
			registerDoParallel(Cl)
		}
		chunkSize  <-  ceiling(n/getDoParWorkers())
		smpopts  <- list(chunkSize = chunkSize)
		
		export = c("bCICore", "rmNA", "rmInf")
	}
	
	
	if (is.vector(X)) { return(bCICore(X, f = f, R = R, conf = conf)) }
	else
	{
		i = NULL
		if (is.data.frame(X))  { X <- factor2matrix(X) }
		
		if (length(which.is.na(X)) > 0)
		{	
			bciObject <- foreach(i = 1:n, .export = export, .options.smp = smpopts, .combine = cbind, .multicombine = TRUE) %dopar%
			bCICore(rmNA(rmInf(X[i,])), f = f, R = R, conf = conf, largeValues = largeValues, skew = skew)
		}
		else		
		{	
			bciObject <- foreach(i = 1:n, .export = export, .options.smp = smpopts, .combine = cbind, .multicombine = TRUE) %dopar%
			bCICore(rmInf(X[i,]),  f = f, R = R, conf = conf, largeValues = largeValues, skew = skew)
		}
		
		stopCluster(Cl)
		bciObject <- t(bciObject)
		colnames(bciObject)[1] = "Estimate"
		return(bciObject)
	}
}


OOBquantiles <- function(object, conf = 0.95, plotting = TRUE, gap = 10, outliersFilter = TRUE, threshold = NULL, thresholdDirection = c("low", "high"))
{
	object = filter.object(object)
	if (is.null(object$forest$OOB.votes)) { stop ("no OOB data. Please enable OOB option when computing randomUniformForest() function") }
	
	X = object$forest$OOB.votes
	B =  ncol(X)
	alpha1 = (1 - conf)/2
	alpha2 = conf + alpha1
	
	Q_1 = apply(X, 1, function(Z) quantile(unique(rmInf(Z)), alpha1))
	Q_2 = apply(X, 1, function(Z) quantile(unique(rmInf(Z)), alpha2))
	SD = apply(X, 1, function(Z) sd(rmInf(Z)))
	predictedObject = data.frame(cbind(object$forest$OOB.predicts, Q_1, Q_2, SD))
	Q1Name = paste("LowerBound", "(q = ", round(alpha1,3), ")", sep ="")
	Q2Name = paste("UpperBound", "(q = ", round(alpha2,3), ")",sep ="")
	colnames(predictedObject) = c("Estimate", Q1Name, Q2Name, "standard deviation")
	
	OOsampleUpper = sum(predictedObject[,3] < object$y)/length(object$y)
	OOsampleLower = sum(predictedObject[,2] > object$y)/length(object$y)
	
	# plot quantiles
	if (plotting == TRUE)
	{
		#require(ggplot2)
		predictedObject2 =  predictedObject
		Y = object$y
		predictedObject2 =  cbind(Y,predictedObject2)
		
		if (!is.null(threshold))
		{ 
			if (thresholdDirection == "low") { 	idx = which(Y <= threshold) }
			else {  idx = which(Y > threshold)  }
			predictedObject2 = predictedObject2[idx, ]
			
			if (length(Y) < 1) { stop("no obervations found using this threshold.") }
		}
		
		if (outliersFilter)
		{
			highOutlierIdx =  which(Y > quantile(Y,0.95))
			lowOutlierIdx =  which(Y < quantile(Y,0.05))
			if (length(highOutlierIdx) > 0 | length(lowOutlierIdx) > 0) 
			{	predictedObject2 = predictedObject2[-c(lowOutlierIdx,highOutlierIdx),]	}
		}		
		predictedObject2 = predictedObject2[seq(1, nrow(predictedObject2), by = gap), ]
		colnames(predictedObject2) = c("Y", "Estimate", "LowerBound", "UpperBound", "standard deviation")
		predictedObject2 = sortMatrix(predictedObject2, 1)
		
		Estimate = UpperBound = LowerBound = NULL
		
		tt <- ggplot(predictedObject2, aes(x = Estimate, y = Y))
			
		plot(tt +  geom_point(size = 3, colour = "lightblue") + geom_errorbar(aes(ymax = UpperBound, ymin = LowerBound)) 
			+ labs(title = "", x = "Estimate (with confidence Interval)", y = "True responses (Y)") )
	}
	
	cat("Probability of being over upper bound (", alpha2*100, "%) of the confidence level: ", round(OOsampleUpper,3) , "\n", sep="")
	cat("Probability of being under the lowest bound (", alpha1*100, "%) of the confidence level: ", round(OOsampleLower,3) , "\n", sep ="")
		
	return(predictedObject)
}

outsideConfIntLevels <- function(predictedObject, Y, conf = 0.95)
{
	alpha1 = (1 - conf)/2
	alpha2 = conf + alpha1
	OOsampleUpper = sum(predictedObject[,3] < Y)/length(Y)
	OOsampleLower = sum(predictedObject[,2] > Y)/length(Y)

	cat("Probability of being over upper bound (", alpha2*100, "%) of the confidence level: ", round(OOsampleUpper,3) , "\n", sep="")
	cat("Probability of being under the lowest bound (", alpha1*100, "%) of the confidence level: ", round(OOsampleLower,3) , "\n", sep ="")
}

scale2AnyValues <- function(X, Y)
{
	Z = vector(length = length(X))
	XX = standardize(X)
	Z = quantile(Y, XX)
	return(Z)
}

perspWithcol <- function(x, y, z, pal, nbColors,...,xlg = TRUE, ylg = TRUE)
{
	# import from http://rtricks.wordpress.com/tag/persp/
	colnames(z) <- y
	rownames(z) <- x

	nrz <- nrow(z)
	ncz <- ncol(z) 

	color <- sort(pal(nbColors), decreasing= TRUE)
	zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
	facetcol <- cut(zfacet, nbColors)
	par(xlog=xlg,ylog=ylg)
	persp(
		as.numeric(rownames(z)),
		as.numeric(colnames(z)),
		as.matrix(z),
		col = color[facetcol],
		...
	)
}

biasVarCov <- function(predictions, target, regression = FALSE, idx = 1:length(target))
{
	if (is.factor(target)) 
	{ 
		target = as.numeric(target) 
		if (length(unique(target)) != 2) { stop("Decomposition currently works only for binary classification.") }
	}
	if (is.factor(predictions)) { predictions = as.numeric(predictions) }
	
	noise = ifelse(length(idx) > 1, var(target[idx]), (target[idx] - mean(target))^2)
	squaredBias = ifelse(length(idx) > 1,(mean(predictions[idx]) - mean(target[idx]))^2,(mean(predictions[idx]) - target[idx])^2) 
	predictionsVar = ifelse(length(idx) > 1, var(predictions[idx]), (predictions[idx] - mean(predictions))^2)
	predictionsTargetCov = ifelse(length(idx) > 1,cov( cbind(predictions[idx], target[idx]) )[1,2], mean( (predictions[idx] - mean(predictions) ) *(target[idx] - mean(target) )))
	
	MSE = noise + squaredBias + predictionsVar - 2*predictionsTargetCov
	
	cat("Noise: ",noise,"\n", sep ="")
	cat("Squared bias: ", squaredBias,"\n", sep ="")
	cat("Variance of estimator: ",  predictionsVar,"\n", sep ="")
	cat("Covariance of estimator and target: ",  predictionsTargetCov, "\n", sep ="")
	
	if (regression)
	{
		object = list(MSE = MSE , squaredBias = squaredBias, predictionsVar = predictionsVar, predictionsTargetCov = predictionsTargetCov)
		cat("\nMSE(Y, Yhat) = Noise(Y) + squaredBias(Y, Yhat) + Var(Yhat) - 2*Cov(Y, Yhat) = ",  MSE, "\n", sep ="")
	}
	else
	{
		object = list(predError = MSE , squaredBias = squaredBias, predictionsVar = predictionsVar, predictionsTargetCov = predictionsTargetCov)
		cat("\nAssuming binary classification with classes {0,1}, where '0' is the majority class.")
		cat("\nMisclassification rate = P(Y = 1)P(Y = 0) + {P(Y = 1) - P(Y_hat = 1)}^2 + P(Y_hat = 0)P(Y_hat = 1) - 2*Cov(Y, Y_hat)")
		cat("\nMisclassification rate = P(Y = 1) + P(Y_hat = 1) - 2*E(Y*Y_hat) = ",  MSE, "\n", sep ="")
	}
	
	return(object)
}
	
dates2numeric <- function(X, whichCol = NULL, inverseDate = FALSE)
{
	lengthChar = nchar(as.character((X[1,whichCol])))

	year = if (inverseDate) { as.numeric(substr(X[,whichCol], 7, 10)) } else { as.numeric(substr(X[,whichCol], 1,4)) }
	month = if (inverseDate) {  as.numeric(substr(X[,whichCol], 4, 5)) } else { as.numeric(substr(X[,whichCol], 6,7)) }
	day =  if (inverseDate) {  as.numeric(substr(X[,whichCol], 1, 2)) } else { as.numeric(substr(X[,whichCol], 9, 10)) }
	
	if (lengthChar > 10) 
	{	
		hour = as.numeric(substr(X[,whichCol], 12,13))	
		if (lengthChar > 13) { minute = as.numeric(substr(X[,whichCol], 15, 16))	}
		else { minute = NULL } 
	}
	else { hour = NULL }
	
	if (is.null(hour)) 	{	return(cbind(year, month, day, X[,-whichCol]))	}
	else {  return(cbind(year, month, day, hour, minute, X[,-whichCol]))	}
	
}	

predictionvsResponses <- function(predictions, responses, plotting = TRUE)
{
	linMod = lm(responses ~ predictions)
	print(summary(linMod))
	if (plotting)
	{
		plot(predictions, responses,  xlab = "Predictions", ylab = "Responses")
		abline(a = linMod$coef[[1]] , b = linMod$coef[[2]], col="green",lwd = 2)
		grid()
	}	
}

rm.tempdir <- function()
{
	path = getwd()
	setwd(tempdir())
	listFiles <- list.files()[grep("file", list.files())]
	if (length(listFiles) > 0) {  file.remove(listFiles) }
	setwd(path)
}

model.stats <- function(predictions, responses, regression = FALSE, OOB = FALSE, plotting = TRUE)
{
    
	if (class(predictions) == "randomUniformForest")
	{
		if (OOB & is.null(predictions$forest$OOB.predicts)) { stop("no OOB predictions in the model.") }
		object = predictions
		if (OOB) { 	predictions = predictions$forest$OOB.predicts }
		else { predictions = predictions$predictionObject$majority.vote }
		
		if (!regression) 
		{
			predictions = as.factor(predictions)
			levels(predictions) = object$classes
		}
	}
	graphics.off()
	#require(pROC)
	if (OOB)  { cat("\nOOB evaluation") }
	else { cat("\nTest set") }
	if (is.factor(predictions) | !regression)
	{
		uniqueResponses = sort(unique(responses))
		if (!is.factor(responses)) { responses = as.factor(responses) }
		confMat = confusion.matrix(predictions, responses)
		#rownames(confMat) = uniqueResponses
		#colnames(confMat)[-ncol(confMat)] = rownames(confMat)
		cat("\nError rate: ")
		testError = generalization.error(confMat)
		cat(round(100*testError, 2), "%\n", sep="")
		cat("\nConfusion matrix:\n")
		print(round(confMat, 4))
		
		if (length(uniqueResponses) == 2)
		{
			AUC = pROC::auc(as.numeric(responses), as.numeric(predictions))
			AUPR = myAUC(as.numeric(predictions), as.numeric(responses), falseDiscoveryRate = TRUE)$auc
			cat("\nArea Under ROC Curve:", round(AUC, 4))
			cat("\nArea Under Precision-Recall Curve:", 
			round(AUPR, 4))
			cat("\nF1-score:", round(fScore(confMat),4))
			if (plotting)
			{
				roc.curve(predictions, responses, uniqueResponses, printValues = FALSE)
				dev.new( )
				roc.curve(predictions, responses, uniqueResponses, falseDiscoveryRate = TRUE, printValues = FALSE)
			}
		}
		cat("\nGeometric mean:", round(gMean(confMat),4),"\n")
		if (nrow(confMat) > 2)
		{ cat("Geometric mean of the precision:", round(gMean(confMat, precision = TRUE),4),"\n") }	
		
		return(list(confusionMatrix = confMat, testError = testError, AUC = if (length(uniqueResponses) == 2) { AUC } else { "none" }, 
		AUPR = if (length(uniqueResponses) == 2) { AUPR } else { "none" }))
	}
	else
	{
		n = length(responses)
		MSE = L2Dist(predictions, responses)/n
		cat("\nMean of squared residuals: ", round(MSE, 6), "\n", sep="") 
		cat("Variance explained: ", round(100*(1 - MSE/var(responses)),2), "%\n\n", sep = "")	
		Residuals <- summary(responses - predictions)
		names(Residuals) = c("Min", "1Q", "Median", "Mean", "3Q", "Max")
		cat("Residuals:\n")
		print(Residuals)
		cat("Mean of absolute residuals: ", round(L1Dist(predictions, responses)/n, 6), "\n", sep="")
		cat("\nLinear modelling:")
		predictionvsResponsesLinMod = predictionvsResponses(predictions, responses, plotting = plotting)
		outputLM = summary(predictionvsResponsesLinMod)
		
		return(list(MSE = MSE, Residuals = Residuals, predictionvsResponsesLinMod = outputLM ))
	}
}

generic.cv <- function(X, Y, nTimes = 1, k = 10, seed = 2014, regression = TRUE, genericAlgo = NULL, specificPredictFunction = NULL, 
metrics = c("none", "AUC", "precision", "F-score", "L1", "geometric mean", "geometric mean (precision)"))
{
	set.seed(seed)
	s = 1
	testError = metricValues = vector()
	n = nrow(X)
	for (i in 1:nTimes)
	{
		# randomize data
		train_test = sample(n, n)
		# divide data in k folds
		foldsIdx = cut(train_test, k, labels = FALSE)
		
		# for each fold (the number 'j'), choose it as a test set 
		# and choose all the others as the training se
		for (j in 1:k)
		{
			# all folds minus the j-th
			trainIdx = which(foldsIdx != j)			
			X1 = X[trainIdx,]
			Y1 = Y[trainIdx]
			
			testIdx = which(foldsIdx == j)
			X2 = X[testIdx,]
			Y2 = Y[testIdx]
									 
			# train and test samples are always converted to an R matrix
			# if classification a factor will be needed
			if (!regression) { Y1 = as.factor(Y1) }
			
			# R fill find genericAlgo() from the global environment
			if (is.null(genericAlgo))
			{ 
				stop("Please provide a wrapper for the algorithm needed to be assessed. For example, if randomForest is needed\n just write : 'genericAlgo <- function(X, Y) randomForest(X, Y)'\n Then, use the option 'genericAlgo = genericAlgo' in generic.cv() function\n")
			}
			else
			{	model = genericAlgo(X1, Y1)	}
			
			# some algorithms like gbm or glmnet require specific arguments
			# for the predict function. One has to create 
			# a specific predict function for them
			if (is.null(specificPredictFunction))
			{ 	predictions = predict(model, X2)	}
			else
			{  predictions = try(specificPredictFunction(model, X2)) }
			
			if (regression)
			{	
				testError[s] = sum((predictions - Y2)^2)/length(Y2)	
			
				if (metrics[1] == "L1") {  metricValues[s] =  sum(abs(predictions - Y2))/length(Y2)	}
			
			}
			else
			{	
				confMat = confusion.matrix(predictions, Y2)
				testError[s] = generalization.error(confMat)
				if (length(unique(Y2)) == 2)
				{	
					if (metrics[1] == "AUC") {  metricValues[s] = pROC::auc(as.numeric(Y2), as.numeric(predictions)) }				
					if (metrics[1] == "precision") 	{ metricValues[s] = confMat[2,2]/sum(confMat[2,]) }				
					if (metrics[1] == "F-score")  	{ metricValues[s] = round(fScore(confMat),4) }
				}
				else 
				{  
					if (metrics[1] == "geometric mean") {   metricValues[s] = round(gMean(confMat),4)	}
					if (metrics[1] == "geometric mean (precision)") {   metricValues[s] = round(gMean(confMat, precision = TRUE),4)	}
				}
			}
			s = s + 1
		}
	}
	
	if (metrics[1] == "none") {  return(list(testError = testError, avgError = mean(testError), stdDev = sd(testError)) ) }
	else {  return(list(testError = testError, avgError = mean(testError), stdDev = sd(testError), metric = metricValues) ) }
}

filterOutliers <- function(X, quantileValue, direction = c("low", "high")) 
if (direction == "low") { X[X < quantileValue] } else { X[X > quantileValue]  }

which.is.nearestCenter <- function(X, Y)
{
	n = nrow(X)
	K = nrow(Y)
	distToY = matrix(NA, n, K)
	for (k in 1:K)	{ distToY[,k] = apply(X, 1, function(Z) L2Dist(Z, Y[k,])) }	
	return(apply(distToY, 1, which.min))
}

variance <- function(X)  (length(X)-1)/length(X)*var(X) 
 
interClassesVariance <- function(X, classes) 
{        
	m <- tapply(X, classes, mean) 
	l <- tapply(X, classes, length)
	return( sum(l* ( m - mean(X) )^2/length(X)) )
} 
 
intraClassesVariance <- function(X, classes) 
{        
	v <- tapply(X, classes, variance) 
	v[is.na(v)] <- 0 
	l <- tapply(X, classes,length) 
	return( sum(l*v/length(X)) )
}
 
mergeOutliers <- function(object)
{
	Z = object$unsupervisedModel$cluster
	if (!is.null(object$unsupervisedModel$clusterOutliers))
	{ Z = c(Z, object$unsupervisedModel$clusterOutliers)	}			
	if (!is.null(names(Z))) 
	{ 
		Z = sortDataframe( data.frame(Z, as.numeric(names(Z))), 2)
		return(Z[,1])
	}
	else
	{	return(Z) }	
}

splitVarCore <- function(X, nsplit, lengthOfEachSplit)
{
	X = as.character(X)
	charLength = nchar(X[1])
	XX = matrix(NA, length(X), nsplit)
	k = 1
	for (j in 1:nsplit)
	{
		for (i in 1:length(X))
		{	XX[i,j] = substr(X[i], k,(k + lengthOfEachSplit[j] - 1)) }
		k = k + lengthOfEachSplit[j] 	
	}

	XX  = fillVariablesNames(XX)
	return(XX)
}

timeStampCore <- function(stamp, n = NULL, begin = 0, end = NULL, windowLength = NULL)
{
	if (is.null(n)) { stop("Please provide 'n', the length of the sample.\n") }
	Z = vector(length = n)
	if (!is.null(end))
	{ Z = seq(begin, end, by = stamp)   }
	else
	{	
		if (!is.null(windowLength))
		{
			K = floor(windowLength/stamp)
			s = 1
			Z[1] = begin
			Z[1:(s + K-1)] = cumsum(c(Z[1], rep(stamp, K-1)))
			s = s + K 
			for (k in 2:(floor(n/K)))
			{
				Z[s:(s + K-1)] = cumsum(rep(stamp, K))			
				s = s + K 
			}
			if (s < n)
			{	
				gapLength = n - s + 1
				Z[s:n] = cumsum(rep(stamp, gapLength))				
			}
		}
		else
		{ stop("Please provide either a 'end' value or a 'windowLength' one.\n") }
	}
	return(Z)	
}

reSMOTE <- function(X, Y, majorityClass = 0, increaseMinorityClassTo = NULL, 
conditionalTo = NULL,
samplingFromGaussian = FALSE,  
samplingMethod = c("uniform univariate sampling", "uniform multivariate sampling", "with bootstrap"),
seed = 2014)
{
	n = nrow(X)
	if (is.data.frame(X)) stop("X must be a matrix.\n")
	if (n != length(Y)) stop("Data and responses don't have same size.\n")
	
	majorityIdx = NULL
	for (i in seq_along(majorityClass))
	{ 	majorityIdx = c(majorityIdx, which(Y == majorityClass[i])) }
	percentMinority = 1 - length(majorityIdx)/n
	minorityIdx = (1:n)[-majorityIdx]
	
	if (is.null(increaseMinorityClassTo)) 
	{ 
		increaseMinorityClassTo = min(2*percentMinority, 0.5)
		cat("Argument 'increaseMinorityClassTo' has not been set.\nOriginal minority class percentage has been raised by a factor 2\n.")
	}
	
	iterations = round(increaseMinorityClassTo/percentMinority, 0)
	
	if (iterations < 1) stop("'increaseMinorityClassTo' must be increased.\n")
	
	if (samplingMethod[1] == "with bootstrap")
	{ 
		cat("'with bootstrap' is only needed as the second argument of 'method' option. Default option will be computed.\n")
		samplingMethod[1] = "uniform univariate sampling"
	}
	
	XXY = vector('list', iterations)
	for (i in 1:iterations)
	{
		XXY[[i]] <- unsupervised2supervised(X, method = samplingMethod[1], conditionalTo = conditionalTo, samplingFromGaussian = samplingFromGaussian, seed = seed,
		bootstrap = if (length(samplingMethod) > 1) { TRUE } else { FALSE })$X[(n+1):(2*n),][minorityIdx,]
	}

	newX = do.call(rbind, XXY)
	rm(XXY)
	newXY = cbind(newX, rep(Y[minorityIdx], iterations))
	p= ncol(newXY)
	colnames(newXY)[p] = "Value"
	if (!is.numeric(Y))
	{	
		levelsY = levels(Y)
		minorityClassIdx = which(levelsY != as.character(majorityClass))
		currentValuesOfY = as.factor(newXY[,p])
		levels(currentValuesOfY) = levelsY[minorityClassIdx]
		newXY = as.data.frame(newXY)
		newXY[,p] = as.factor(currentValuesOfY)
		colnames(newXY)[p] = "Class"
		X = as.data.frame(X)
		XY = (cbind(X,Y))
		colnames(XY)[p] = "Class"
		ZY  = rbind(XY, newXY)
	}
	else
	{
		XY = (cbind(X,Y))
		colnames(XY)[p] = "Value"
		ZY  = rbind(XY, newXY)	
	}
	cat("Proportion of examples of minority class in the original sample: ", round(100*(length(minorityIdx))/n, 2), "%.\n", sep="")
	cat("Proportion of examples of minority class in the new sample: ", round(100*length(which(ZY[,"Class"] != majorityClass))/nrow(ZY), 2), "%.\n", sep="")
	cat(nrow(ZY), "instances in the new sample.\n")
	
	return(ZY)
} 
# END OF FILE