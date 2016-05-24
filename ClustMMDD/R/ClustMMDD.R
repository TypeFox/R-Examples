#######################################
####	Main R functions interfacing R and cpp 
#######################################
## library("Rcpp")
## sourceCpp("ClustMMDD.cpp")

#########################################
###	Other useful functions
#########################################

#########################################
##	Transform a data.frame from R to Rcpp
#########################################

#'	Transform a data.frame for \code{\link{ClustMMDD}} main functions compatibility.
#'	@title	Transform a data.frame for \code{\link{ClustMMDD}} main functions compatibility
#'
#'	@description
#'	\code{dataR2C(x, ploidy)} returns a list. 
#'
#'	@param	x	A data.frame or a matrix with number of columns equal to \code{ncol(x)*ploidy}.
#'	@param	ploidy	The number of columns for each variable. For example, ploidy = 2 for genotypic data from diploid individual.

#'	@return	A list of elements needed for \code{\link{ClustMMDD}} main functions :
#'	\itemize{
#'		\item{$data}{A matrix compatible with \code{\link{ClustMMDD}} main functions.}
#'		\item{$ploidy}{The number of columns for each variable.}
#'		\item{$N}{The number of lines in x.}
#'		\item{$P}{The number of variables equal to ncol(x)/ploidy.}
#'		\item{$N_LEVELS}{The vector of the numbers of levels for the variables.}
#'		\item{$LEVELS}{The levels for the variables.}
#'		\item{$COUNT}{The observed counts of the levels.}
#'		\item{$FREQ}{The observed frequencies.}
#'	}
#'	
#'	@author	Wilson Toussile
#'
#'	@seealso \code{\link{cutEachCol}} in \code{\link{ClustMMDD}} package.
#'
#'	@examples
#'	\dontrun{
#'	data(genotype1)
#'	head(genotype1) # Each data is a set of 2 observed alleles. For example "101102" = ("101", "102").
#'	genotype2 = cutEachCol(genotype1[, -11], 2) #ploidy = 2 for these genotype data
#'	head(genotype2)
#'	genotype3 = dataR2C(genotype2, ploidy = 2)
#'	head (genotype3$data)
#'	str(genotype3)
#'	}
#'	@export
dataR2C <- function(x, ploidy = 1)
{
  if(!is.element(ploidy, 1:5)) stop("ploidy not in 1:5")
  if(ncol(x)%%ploidy != 0) stop("Number of columns in data not corresponding to ploidy !")

  x = as.matrix(x)
  N = nrow(x)	
  P = ncol(x)/ploidy 
  
  if(ploidy == 1) 
  {
	x1 = data.frame(x)
  }
  else
  {
	x1 = sapply(1:P, function(j)
	{
	  j0 = (j-1)*ploidy + 1
	  return(x[, j0:(j0+ploidy-1)])
	})
  }
  
  COUNT = lapply(1:P, function(j) return(table(unlist(x1[, j]))))
  
  xC = matrix(match(x1[ ,1], names(COUNT[[1]])), ncol = ploidy)
  if(P>1)
  {
	for(j in 2:P)
	  xC = cbind(xC, matrix(match(x1[ ,j], names(COUNT[[j]])), ncol = ploidy))
  }
  
  Freq = lapply(1:P, function(j) return(COUNT[[j]]/(ploidy*N)))
  
  N_levels = sapply(1:P, function(j) return(length(COUNT[[j]])))
  
  
  LevelsInR = lapply(1:P, function(j) return(as.character(names(COUNT[[j]]))))
  LevelsR2C = unlist(LevelsInR)

  xC1 = xC - 1
  
  return(list(data = as.matrix(xC1), 
			  ploidy = as.integer(ploidy), 
			  N = as.integer(N), 
			  P = as.integer(P),
			  N_LEVELS = as.integer(N_levels), 
			  LEVELS = as.character(LevelsR2C), 
			  LEVELS_R = LevelsInR, 
			  COUNT = as.integer(unlist(COUNT)), 
			  FREQ = as.double(unlist(Freq))))
}

## cutEachCol
#'	Cut each column of a data.frame or matrix of strings to a given number of columns.
#' @description
#'	Each string in the data.frame represents a set of \code{\link{ploidy}} observations from the same set of levels. For example "101102" may represents {"101", "102"}, in which case ploidy = 2.
#'	@title	cutEachCol
#'
#'	@param	xdata	A data.frame or a matrix of strings.
#'	@param	ploidy	The number of columns for each column. For example, for genotypic data from diploid individual, \code{ploidy = 2} : for example, a data "ab" represents \code{\{"a", "b"\}} observed alleles.
#'
#'	@return	A list of elements needed for \code{\link{ClustMMDD}} main functions :
#'	\itemize{
#'		\item{data}{A matrix compatible with \code{\link{ClustMMDD}}  main functions.}
#'		\item{ploidy}{The number of columns for each column.}
#'		\item{N}{The number of lines in xdata.}
#'		\item{P}{The number of columns in xdata}
#'		\item{N_LEVELS}{The vector of the numbers of levels.}
#'		\item{LEVELS}{The levels for the variables.}
#'		\item{COUNT}{The observed counts of the levels.}
#'		\item{FREQ}{The observed frequencies.}
#'	}
#'
#'
#'	@author	Wilson Toussile
#'
#'	@seealso \code{\link{dataR2C}} for \code{\link{ClustMMDD}} data format.
#'
#'	@examples
#'	\dontrun{
#'	data(genotype1)
#'	head(genotype1) # Each data is a set of 2 observed alleles. For example "101102" = ("101", "102").
#'	genotype2 = cutEachCol(genotype1[,-11], ploidy = 2) #ploidy = 2 for these genotype data
#'	head(genotype2)
#'	}
#'	@export
cutEachCol <- function(xdata, ploidy)
{
  if(!is.data.frame(xdata))
	stop("Argument xdata must be a data.frame")
	
  NberMax = 5
  if(!is.element(ploidy, 1:NberMax)) 
		stop(paste("ploidy not in ", 1:NberMax))
		
  n = nrow(xdata)
  
  xdata = matrix(as.character(unlist(xdata)), nrow = n)
  
	##Rcpp
  return(try(cutEachColInN(xdata, ploidy)))
}

## Set EM options
#'	Set Expectation and Maximization options.
#'
#'	@title	Set Expectation and Maximization options.
#'
#'	@param	emOptions	A list of options needed by the Expectation and maximization algorithm :
#'	\itemize{
#'		\item{epsi}{In [1e-5, 1e-20], it is upper bound of the relative increasing on log-likelihood.} 
#'		\item{typeSmallEM}{In c(0, 1, 2) : 0 = classic EM, 1 = SEM, 2 = CEM.}
#'		\item{typeEM}{In c(0, 1, 2) : 0 = classic EM, 1 = SEM, 2 = CEM.} 
#'		\item{nberSmallEM}{The number of random parameter point from which to perform a certain \code{nberIterations} of EM runs.} 
#'		\item{nberIterations}{The number of iterations for each small EM.} 
#'		\item{nberMaxIterations}{The maximum number of iterations if EM algorithm converge hardly.} 
#'		\item{putThreshold}{If \code{TRUE}, the probabilities of levels are assumed to be positive in all clusters.}
#'	}
#'
#'	@details
#'	Use \code{setEmOptions()} to set all options to default.
#'	
#'	@author	Wilson Toussile.
#'
#'	@references
#'	\itemize{
#'		\item  \href{http://projecteuclid.org/euclid.ejs/1379596773}{Dominique Bontemps and Wilson Toussile (2013)} : Clustering and #'	variable selection for categorical multivariate data. Electronic Journal of Statistics, Volume 7, 2344-2371, ISSN.
#'		\item 	\href{http://www.springerlink.com/content/r7k02q167120xl64/}{Wilson Toussile and Elisabeth Gassiat (2009)} : Variable #'	selection in model-based clustering using multilocus genotype data. Adv Data Anal Classif, Vol 3, number 2, 109-134.
#'	}
#'
#'	@seealso	\code{\link{EmOptions}} for getting the current EM options.
#'
#'	@keywords	EM
#'
#'	@examples
#'	\dontrun{
#'	EmOptions()
#'	setEmOptions(list(epsi = 1e-6))
#'	EmOptions()
#'	setEmOptions()	# To set default values
#'	EmOptions()
#'	}
#'	@export
setEmOptions <- function(emOptions = list(epsi = NULL, 
										typeSmallEM = NULL, 
										typeEM = NULL, 
										nberSmallEM = NULL, 
										nberIterations = NULL, 
										nberMaxIterations = NULL, 
										putThreshold = NULL))
{
	##Rcpp
	try(
  	EmSettings(xepsi = ifelse(is.null(emOptions[["epsi"]]), -1.0, emOptions[["epsi"]]),
  	xnberSmallEM = ifelse(is.null(emOptions[["nberSmallEM"]]), -1, emOptions[["nberSmallEM"]]),
  	xnberIterations = ifelse(is.null(emOptions[["nberIterations"]]), -1, emOptions[["nberIterations"]]),
  	xtypeSmallEM = ifelse(is.null(emOptions[["typeSmallEM"]]), 0, emOptions[["typeSmallEM"]]),
  	xtypeEM = ifelse(is.null(emOptions[["typeEM"]]), 0, emOptions[["typeEM"]]),
  	xnberIterLongEM = ifelse(is.null(emOptions[["nberMaxIterations"]]), -1, emOptions[["nberMaxIterations"]]),
  	xputThreshold = ifelse(is.null(emOptions[["putThreshold"]]), FALSE, emOptions[["putThreshold"]])))
}

##	Get the current EM options
#'	@title	Display current Expectation and Maximization options.
#'
#'	@description
#'	Display the Expectation and Maximization algorithm current options.
#'
#'	@return	A list of EM options :
#'	\itemize{
#'		\item{epsi : }{The upper bound of the relative increasing on log-likelihood.} 
#'		\item{nberSmallEM : }{The number of random parameters points from which to run small EM. The estimated parameter point associated to the 
#'		maximum log-likelihood is then used for a final EM run.}
#'		\item{nberIterations : }{The number of iterations in each small EM.}
#'		\item{typeSmallEM : }{0 = classic EM, 1 = SEM and 2 = CEM.}
#'		\item{typeEM : }{0 = classic EM, 1 = SEM and 2 = CEM.}
#'		\item{nberMaxIterations : }{The maximum number of iterations in the final EM if the convergence is slow.}
#'		\item{putThreshold : }{The indication whether all parameter estimates are positive.}
#'	}
#'
#'	@author	Wilson Toussile.
#'
#'	@references
#'	\itemize{
#'		\item  \href{http://projecteuclid.org/euclid.ejs/1379596773}{Dominique Bontemps and Wilson Toussile (2013)} : Clustering and #'	variable selection for categorical multivariate data. Electronic Journal of Statistics, Volume 7, 2344-2371, ISSN.
#'		\item 	\href{http://www.springerlink.com/content/r7k02q167120xl64/}{Wilson Toussile and Elisabeth Gassiat (2009)} : Variable #'	selection in model-based clustering using multilocus genotype data. Adv Data Anal Classif, Vol 3, number 2, 109-134.
#'	}
#'
#'	@keywords EM SEM CEM
#'
#'	@seealso \code{\link{setEmOptions}} for set EM options.
#'	
#'	@examples
#'	\dontrun{
#'	EmOptions()
#'	setEmOptions(list(epsi = 1e-6))
#'	EmOptions()
#'	setEmOptions()	# To set default values
#'	EmOptions()
#'	}
#'	@export
EmOptions <- function()
{
	##Rcpp
	return (getEmOptions_Rcpp())
}

##	A general function for the EM algorithm ##
#'	@title	Perform the Expectation and Maximization algorithm followed by the maximum a posteriori clustering.
#'
#'	@description
#'	Compute an approximation of the maximum likelihood parameter estimates using Expectation and Maximization (EM) algorithm.
#'
#'	@param	xdata	A matrix of strings with the number of columns equal to ploidy * (number of variables).
#'	@param	K	The number of clusters (or populations).
#'	@param	S	The subset of clustering variables on the form of vector of logicals indicating the selected variables of clustering.
#'	@param	emOptions	A list of EM options (see \code{\link{EmOptions}} and \code{\link{setEmOptions}}).
#'	@param	ploidy	The number of occurrences for each variable in the data. For example, ploidy = 2 for genotype data from diploid individual.
#'	@param	cte	A double used for the selection criterion named \code{CteDim} in which the penalty function is \code{pen(K,S)=cte*dim}, 
#'	where \code{dim} is the number of free parameters.
#'
#'	@return	A list of 
#'	\itemize{
#'		\item{$N}{The size (number of lines) of the dataset.}
#'		\item{K}{The number of clusters (populations).}
#'		\item{S}{A vector of logicals indicating the selected variables for clustering.}
#'		\item{dim}{The number of free parameters.}
#'		\item{pi_K}{The vector of mixing proportions.}
#'		\item{prob}{A list of matrices, each matrix being the probabilities of a variables.}
#'		\item{logLik}{The log-likelihood.}
#'		\item{entropy}{The entropy.}
#'		\item{criteria}{Criteria values c(BIC, AIC, ICL, CteDim).}
#'		\item{Tik}{A stochastic matrix given belonging probabilities of individuals in the sample.}
#'		\item{mapClassif}{Maximum a posteriori classification.}
#'		\item{NbersLevels}{The numbers of observed levels.}
#'		\item{levels}{The observed levels.}
#'	}
#'	
#'
#'	@author	Wilson Toussile.
#'
#'	@references
#'	\itemize{
#'		\item  \href{http://projecteuclid.org/euclid.ejs/1379596773}{Dominique Bontemps and Wilson Toussile (2013)} : Clustering and #'	variable selection for categorical multivariate data. Electronic Journal of Statistics, Volume 7, 2344-2371, ISSN.
#'		\item 	\href{http://www.springerlink.com/content/r7k02q167120xl64/}{Wilson Toussile and Elisabeth Gassiat (2009)} : Variable #'	selection in model-based clustering using multilocus genotype data. Adv Data Anal Classif, Vol 3, number 2, 109-134.
#'	}
#'
#'	@keywords	EM Classification Clustering
#'
#'	@seealso \code{\link{dataR2C}}.
#'
#'	@examples
#'	\dontrun{
#'	data(genotype1)
#'	head(genotype1) # Each data is a set of 2 observed alleles. For example "101102" = ("101", "102").
#'	genotype2 = cutEachCol(genotype1[, -11], ploidy = 2) #ploidy = 2 for these genotype data
#'	head(genotype2)
#'	
#'	#See the EM options
#'	EmOptions() # Options can be set by \code{\link{setEmOptions()}}
#'	par5 = em.cluster.R (genotype2, K = 5, S = c(rep(TRUE, 8), rep(FALSE, 2)), ploidy = 2)
#'	str(par5)
#'	head(par5$Tik)
#'	}
#'	@export
em.cluster.R <- function(xdata, K, S,
				 ploidy = 1,
				 emOptions = list(
				 epsi = NULL,
				 typeSmallEM = NULL,
				 typeEM = NULL, 
				 nberSmallEM = NULL,
				 nberIterations = NULL,
				 nberMaxIterations = NULL,
				 putThreshold = NULL),
				 cte = 1.0)
{
  K = as.integer(K)
  S = as.logical(S)
  ploidy = as.integer(ploidy)
  
  ##	Transform data for Rcpp function
  data2C = dataR2C(xdata, ploidy = ploidy)
  
  n = as.integer(data2C$N)
  p = as.integer(ncol(xdata)/ploidy)
  
  ##	Compatibility of dataset dimensions
  if(p != length(S)) stop("Incompatibility between the number of variables and the length of S !")
  
  S = as.logical(S)
  n_levels <- data2C$N_LEVELS
  count <- data2C$COUNT
  levels <- data2C$LEVELS
  freq <- data2C$FREQ
  
  setEmOptions(emOptions)
  
  outPar = try(EM_Rcpp(t(data2C$data), ploidy , levels, n_levels, count, freq, K, S, cte))
    
  try(initialiseEmSettings())
  
  invisible(setModelKS(outPar))  
}

##	Selection of the number of components for a given S
#'	Perform a selection of the number K of clusters for a given subset S of clustering variables.
#'
#'	@title	Selection of the number K of clusters.
#'
#'	\code{selectK.R(xdata, S, Kmax)} returns a list of estimated paramaters for each selection criteria.
#'
#'	@param	xdata	A dataset.
#'	@param	S	A subset of clustering variables on the form of logical vector of the same length P as the number of variables in x.
#'	@param	Kmax	The maximum number of clusters to be explored.
#'	@param	Kmin	The minimum number of clusters to be explored. The default value is set to 1.
#'	@param	emOptions	A list of EM options (see \code{\link{EmOptions}} and \code{\link{setEmOptions}}).
#'	@param	ploidy	The number of occurrences for each variable in the data. For example, ploidy = 2 for genotype 
#'	@param	cte	A double used for the selection criterion named \code{CteDim} in which the penalty function is 
#'	\code{pen(K,S)=cte*dim}, where \code{dim} is the number of free parameters.
#'
#'	@return	A list of estimated paramaters for each selection criteria. 
#'
#'	@keywords	EM Classification	Clustering
#'
#'	@author	Wilson Toussile
#'
#'	@references
#'	\itemize{
#'		\item  \href{http://projecteuclid.org/euclid.ejs/1379596773}{Dominique Bontemps and Wilson Toussile (2013)} : Clustering and #'	variable selection for categorical multivariate data. Electronic Journal of Statistics, Volume 7, 2344-2371, ISSN.
#'		\item 	\href{http://www.springerlink.com/content/r7k02q167120xl64/}{Wilson Toussile and Elisabeth Gassiat (2009)} : Variable #'	selection in model-based clustering using multilocus genotype data. Adv Data Anal Classif, Vol 3, number 2, 109-134.
#'	}
#'
#'	@seealso \code{\link{backward.explorer}} for more exploration of the competing models, \code{\link{dimJump.R}} 
#'	for data drven calibration of the penality function, and \code{\link{selection.model.R}} for model selection.
#'
#'	@examples
#'	\dontrun{
#'	data(genotype1)
#'	head(genotype1) # Each data is a set of 2 observed alleles. For example "101102" = ("101", "102").
#'	genotype2 = cutEachCol(genotype1[, -11], ploidy = 2) #ploidy = 2 for these genotype data
#'	head(genotype2)
#'	S = c(rep(TRUE, 8), rep(FALSE, 2))
#'
#'	outPut = selectK.R(genotype2, S, Kmax = 6, ploidy = 2, Kmin=1)
#'	outPut[["BIC"]]
#'	}
#'	@export
selectK.R <- function(xdata, S, Kmax,
					  ploidy = 1, 
					  Kmin=1, 
					  emOptions = list(
					  epsi = 1e-5,
					  nberSmallEM = 20,
					  nberIterations = 15,
					  nberMaxIterations = 5000,
					  typeSmallEM = 0,
					  typeEM = 0,
					  putThreshold = FALSE),
					  cte = 1.0,
					  project = deparse(substitute(xdata)))
{
  p = as.integer(ncol(xdata)/ploidy)
  if(length(S) != p) stop("Incorrect dimensions !")
  if(Kmax<=Kmin) stop("Incorrect range of the number of components: choose Kmax >= Kmin !")
  
  data2C = dataR2C(xdata, ploidy = ploidy)
  ploidy = as.integer(ploidy)
  S = as.logical(S)
  n = data2C$N
  p = data2C$P
  count <- data2C$COUNT
  levels <- data2C$LEVELS
  n_levels <- data2C$N_LEVELS
  freq <- data2C$FREQ
  
  setEmOptions(emOptions)
	NBER_CRITERIA = getNberCriteria_Rcpp()
  
  criteria = matrix(ncol = NBER_CRITERIA)
  outPar = list()
  
  FileNameExploredModels = paste(project, "_", "ExploredModels.txt", sep="")
	if(!file.exists(FileNameExploredModels))
	{
  	file.create(FileNameExploredModels)
		cat(c("N", "P",  "K", paste("S", 1:p, sep=""), "logLik",  "dim", "entropy\n"), file = FileNameExploredModels)
	}
  
  #FIXME TFaire en sorte que seuls qui ne sont pas dans le fichier le soient.
  for(K in Kmin:Kmax)
  {
    cat("\n K = ", K ,"\n")

		K = as.integer(K)
	
		##Rcpp
		outParK = try(EM_Rcpp(t(data2C$data), ploidy , levels, n_levels, count, freq, K, S, cte))

		outPar[[K - Kmin + 1]] = outParK
		criteria = rbind(criteria, outParK[["criteria"]])

		##Rcpp
		xfound = isInFile_Rcpp(K, outParK[["S"]], FileNameExploredModels, header = TRUE) #list(TrueFalse, line, dim, logLik)

		if(!xfound[["TrueFalse"]])
		{
		  outParKx = list(N = n, P = p, K = K, S = as.logical(outParK[["S"]]),
		  logLik = as.double(outParK[["logLik"]]), dim = as.integer(outParK[["dim"]]),
		  entropy = as.double(outParK[["entropy"]]))

		  writeModelInFile_Rcpp(outParKx, FileNameExploredModels)
		}
	}
  
  criteria = criteria[-1, ]
  colnames(criteria) = getCriteriaNames_Rcpp()
  
  #Best model by criterion
  bestK = apply(criteria, 2, which.min)
  outPutList = list()
  for(i in 1:NBER_CRITERIA)
		outPutList[[i]] = setModelKS(outPar[[bestK[i]]])
	
  names(outPutList) = getCriteriaNames_Rcpp()
  
  setEmOptions()
  
  invisible(outPutList)
}


##	Find a model in a file
#'	Find a given model defined by (K, S) in a file.
#'	@title	Find a model in a file.
#'
#'	\code{isInFile.R(K, S, file, header)}
#'	@param	K	A positive integer which is the number of clusters.
#'	@param	S	A vector of logicals of the length equal to the number of variables indicating the clustering variables. 
#'	@param	file	A file where to find the model.
#'	@param	header	A logical indicating if the file containts header or not.
#'
#'	@return	0 A list : 
#'	\itemize{
#'		\item A logical indicating if the given model was found and the following if TRUE.
#'		\item The line where the given is in the file.
#'		\item The size of the dataset from which the model was estimated.
#'		\item The log-likelihood.
#'		\item The dimension of the model.
#'		\item The entropy associated to estimated parameters of the models.
#'	}
#'
#'	@author	Wilson Toussile
#'
#'	@references
#'	\itemize{
#'		\item  \href{http://projecteuclid.org/euclid.ejs/1379596773}{Dominique Bontemps and Wilson Toussile (2013)} : Clustering and #'	variable selection for categorical multivariate data. Electronic Journal of Statistics, Volume 7, 2344-2371, ISSN.
#'		\item 	\href{http://www.springerlink.com/content/r7k02q167120xl64/}{Wilson Toussile and Elisabeth Gassiat (2009)} : Variable #'	selection in model-based clustering using multilocus genotype data. Adv Data Anal Classif, Vol 3, number 2, 109-134.
#'	}
#'	
#'	@examples
#'	\dontrun{
#'	data(genotype1)
#'	head(genotype1) # Each data is a set of 2 observed alleles. For example "101102" = ("101", "102").
#'	genotype2 = cutEachCol(genotype1[, -11], ploidy = 2) #ploidy = 2 for these genotype data
#'	head(genotype2)
#'	S = c(rep(TRUE, 8), rep(FALSE, 2))
#'
#'	outPut = selectK.R(genotype2, S, Kmax = 6, ploidy = 2, Kmin=1)
#'	isInFile.R(K = 5, S, "genotype2_ExploredModels.txt", header = TRUE)
#'	isInFile(K = 5, rep(TRUE, 10), "genotype2_ExploredModels.txt", header = TRUE)
#'	}
#'	@export
isInFile.R <- function(K, S, file, header)
{
  ## isInFile_Rcpp() from Cpp returns list(TrueFalse, line, dim, logLik)
  if(!is.logical(S))
		stop("Argument S must be a vector of logicals")
		
	outList = isInFile_Rcpp(K, S, file, header)
	outList[["line"]] = outList[["line"]] + 1
  return(outList)
}


##	Read or compute the EM 
##	TODO : Write a cpp version of this function in the future
##	Not available for users
#'	Read a given model from a file or compute if not found. This function is not available for users.
#'	@title	Read a given model from a file or compute the estimates if not found.
#'
#'	@param	xdata	A list of dataset and several description paramaters seach as frequencies.
#'	@param	xK	The number of components (clusters or populations).
#'	@param	xS	The subset of relevant variables.
#'	@param	xReferenceModel	The indicator of if the model is a reference model in an exclusion step of the backward-stepwise explorer.
#'	@param	xReferenceModelsIndex	The vector indicating the models that have been at once reference at an exclusion step.
#'	@param	xNberExploredModels	The current number of explored models.
#'	@param	xFileNameExploredModels	The explored models.
#'	@param	cte	A constant real.
#'	@param	header	Indication of the presence of header in the file.
#'
#'	@author	Wilson Toussile
#'
#'	@references
#'	\itemize{
#'		\item  \href{http://projecteuclid.org/euclid.ejs/1379596773}{Dominique Bontemps and Wilson Toussile (2013)} : Clustering and #'	variable selection for categorical multivariate data. Electronic Journal of Statistics, Volume 7, 2344-2371, ISSN.
#'		\item 	\href{http://www.springerlink.com/content/r7k02q167120xl64/}{Wilson Toussile and Elisabeth Gassiat (2009)} : Variable #'	selection in model-based clustering using multilocus genotype data. Adv Data Anal Classif, Vol 3, number 2, 109-134.
#'	}
#'
#'	@seealso \code{\link{dimJump.R}} for data drven calibration of the penality function, and 
#'	\code{\link{selection.model.R}} for model selection.
#'
read.or.compute <- function(xdata, 
							xK, 
							xS, 
							xReferenceModel, 
							xReferenceModelsIndex,
							xNberExploredModels,
							xFileNameExploredModels,
							cte = as.double(1.0),
							header = TRUE)
{
  LocalxReferenceModelsIndex <- deparse(substitute(xReferenceModelsIndex))
  LocalNberExploredModels <- deparse(substitute(xNberExploredModels))
  
  ##	Find the model in xFileNameExploredModels
  xfound = isInFile_Rcpp(xK, xS, xFileNameExploredModels, header) #list(TrueFalse, line - 1, dim, logLik) 
  
  ## debug
  #print(xfound)
  
  xN = as.integer(xdata$N)
	NBER_CRITERIA = getNberCriteria_Rcpp()
  criteria = double(NBER_CRITERIA)
  
  if(xfound[["TrueFalse"]])
  {
		numModel = xfound[["line"]] + 1
		if(xReferenceModel) xReferenceModelsIndex[numModel] = TRUE
	  
		##Rcpp::DoubleVector computeCriteria_Rcpp(double lv, int dim, int N, double entropy = 0, double cte = 1.0) from cpp.
		criteria = computeCriteria_Rcpp(xfound[["logLik"]], xfound[["dim"]], xN, xfound[["entropy"]], cte)
		return(criteria)
  }
  else
  {
		xP <- as.integer(xdata$P)
		xploidy <- as.integer(xdata$ploidy)
		xN_LEVELS = as.integer(xdata$N_LEVELS)
		xLEVELS <- as.integer(xdata$LEVELS)
		xCOUNT <- as.integer(xdata$COUNT)
		xFREQ <- as.double(xdata$FREQ)
		xK <- as.integer(xK)
		xS <- as.logical(xS)
	
		##Rcpp::List EM(Rcpp::IntegerMatrix tab, int n_occurrences, Rcpp::CharacterVector levels, Rcpp::IntegerVector n_levels, Rcpp::IntegerVector levels_count,
		#Rcpp::DoubleVector levels_freq, int K, Rcpp::LogicalVector S, double cte = 1.0)
		outPar <- try(EM_Rcpp(t(xdata$data), xploidy , xLEVELS, xN_LEVELS, xCOUNT, xFREQ, xK, xS, cte))
	
		## Now, the number of models = number of lines in FileNameExploredModels
		xNberExploredModels <- nberOfLines(xFileNameExploredModels) + 1
		xReferenceModelsIndex[xNberExploredModels] <- xReferenceModel
		
		outParX = list(N = xN, P = xP, K = xK, S = xS, logLik = outPar[["logLik"]], dim = outPar[["dim"]], entropy = outPar[["entropy"]])
	
		##Rcpp
		writeModelInFile_Rcpp(outParX, xFileNameExploredModels)
	
		assign(LocalxReferenceModelsIndex, xReferenceModelsIndex, pos = parent.frame())
		assign(LocalNberExploredModels, xNberExploredModels, pos = parent.frame())
	
		return(outPar[["criteria"]])
  }  
}

##	Backward explorer
#'	Gather a set of most competive models using a backward-stepwise strategy.
#'	@title	Gather a set most competive models.
#'
#'	@param	x	The data matrix.
#'	@param	Kmax	The maximum number of clusters to be explored.
#'	@param	Criterion	The model selection criterion in c("BIC", "AIC", "ICL", "CteDim") used for exploration (see details).
#'	@param	ForceExclusion	The indication of whether to force exclusion or not. The default value is set to FALSE.
#'	@param	Kmin	The minimum number of clusters. The default value is set to 1.
#'	@param	Smin	A logical vector that indicates variables to include in the set of selected variables. 
#'	The default value NULL: no variable is preselected.
#'	@param	ploidy	The number of columns for each variable in the data. For example, ploidy = 2 for genotypic data from diploid individual.
#'	@param	emOptions	A list of EM options (see \code{\link{EmOptions}} and \code{\link{setEmOptions}}).
#'	@param	project	The name of the project. The default value is the name of the dataset.
#'
#'	@return A data.frame of selected models for the proposed criteria.
#'	
#'	@details 1
#'	If the penalized log-likelihood criteria is \code{CteDim}, a sequence of penalty functions of the form
#'	\code{CteDim*dim} is used. In this shape of penalty function, \code{CteDim} is in \code{[0.5, log(N)]}, where 
#'	\code{N} is the number of individuals in the sample data. Thus, \code{AIC} and \code{BIC} penalties are i the sequence
#'	of penalties.
#'  
#'	@author	Wilson Toussile
#'
#'	@references
#'	\itemize{
#'		\item  \href{http://projecteuclid.org/euclid.ejs/1379596773}{Dominique Bontemps and Wilson Toussile (2013)} : Clustering and #'	variable selection for categorical multivariate data. Electronic Journal of Statistics, Volume 7, 2344-2371, ISSN.
#'		\item 	\href{http://www.springerlink.com/content/r7k02q167120xl64/}{Wilson Toussile and Elisabeth Gassiat (2009)} : Variable #'	selection in model-based clustering using multilocus genotype data. Adv Data Anal Classif, Vol 3, number 2, 109-134.
#'	}
#'	
#'	@keywords	Backward-stepwise	Cluster	Classification
#'
#'	@seealso \code{\link{dimJump.R}} for data drven calibration of the penality function, and 
#'	\code{\link{selection.model.R}} for model selection.
#'
#'	@examples
#'	\dontrun{
#'	data(genotype1)
#'	head(genotype1) # Each data is a set of 2 observed alleles. For example "101102" = ("101", "102").
#'	genotype2 = cutEachCol(genotype1[, -11], ploidy = 2) #ploidy = 2 for these genotype data
#'	########
#'	# The following command create a file "genotype2_ExploredModels.txt" that containts the explored models.
#'	
#'	#output = backward.explorer(genotype2, Kmax = 10, ploidy = 2, Kmin = 1, Criterion = "CteDim")
#'	
#'	
#'	data(genotype2_ExploredModels)
#'	head(genotype2_ExploredModels)
#'	}
#'	@export
backward.explorer <- function(x,
									 Kmax, 
									 Criterion,
									 ploidy = 1,
									 ForceExclusion = FALSE,
									 emOptions = list(
									 epsi = NULL,
									 nberSmallEM = NULL,
									 nberIterations = NULL,
									 nberMaxIterations = NULL,
									 typeSmallEM = NULL,
									 typeEM = NULL,
									 putThreshold = NULL),
									 Kmin = 1, 
									 Smin = NULL,
									 project = deparse(substitute(x)))
{  
  # Kmax > Kmin
  if(Kmax == 1 | Kmax < Kmin) 
	stop("Incorrect range of the number of components: choose Kamx > 1 and Kmax >= Kmin !")
	
  ##	Number of criteria
  ##Rcpp
	NBER_CRITERIA = getNberCriteria_Rcpp()
  NAMES_CRITERIA = getCriteriaNames_Rcpp()
	
  if(!is.element(Criterion, NAMES_CRITERIA)) 
    stop("Choose the argument Criterion in c(BIC, AIC, ICL, CteDim) !")

  CRIT = match(Criterion, NAMES_CRITERIA)

  ##Set grid for exploration
  if(CRIT == 4) cte = as.double(seq(from = 0.5, to = log(1000), by = 0.5))
  if(CRIT != 4)   cte = as.double(1.0)
	
  ##	Output files
  #	For explored models
  FileNameExploredModels = paste(project, "_", "ExploredModels.txt", sep="")

  ##	Transform data
  data2C = dataR2C(x, ploidy = ploidy)
  
  P = data2C$P
  
  if(!file.exists(FileNameExploredModels))
	{
  	file.create(FileNameExploredModels)
		cat(c("N", "P",  "K", paste("S", 1:P, sep=""), "logLik",  "dim", "entropy\n"), file = FileNameExploredModels)
	}
  
  ##	Consider all variables if is.null(Smin)
  if(is.null(Smin)) Smin = rep(FALSE, P)
  else if(all(Smin)) stop("Incorrect Smin !")
  
  
  ##	
  ploidy = as.integer(ploidy)
  xdata = data2C$data
  N = data2C$N
  P = data2C$P
  COUNT <- data2C$COUNT
  LEVELS <- data2C$LEVELS
  N_LEVELS <- data2C$N_LEVELS
  FREQ <- data2C$FREQ
  
  ##	Explored models on matrix form 
  NberExploredModels = 0
  ReferenceModelsIndex = c()
  
  ##	If FileNameExploredModels is not empty, initialise ReferenceModelsIndex FIXME
  if(file.exists(FileNameExploredModels))
  {
		NberExploredModels = nberOfLines(FileNameExploredModels)
		ReferenceModelsIndex = rep(FALSE, NberExploredModels)
  }
  
  ##	Control parameters
  HasExcluded = FALSE
  HasIncluded = FALSE
    
  if(Kmin == 1)
  {
		S = rep(FALSE, P)
	
		xfound = isInFile_Rcpp(Kmin, S, FileNameExploredModels, header = TRUE) #list(TrueFalse, line, dim, logLik)
		if(!xfound[["TrueFalse"]])
		{
			##Rcpp
	  	outParK = try(EM_Rcpp(t(xdata), ploidy , LEVELS, N_LEVELS, COUNT, FREQ, Kmin, S, 1.0))
			
	  	outParKX = list(N = N, P = P, K = Kmin, S = S, logLik = outParK[["logLik"]], dim = outParK[["dim"]], entropy = outParK[["entropy"]])
	  	writeModelInFile_Rcpp(outParKX, FileNameExploredModels)
	  	NberExploredModels = NberExploredModels + 1
	  	ReferenceModelsIndex[NberExploredModels] <- FALSE
		}
		else
		{
	  	ReferenceModelsIndex[xfound[["line"]] + 1] <- FALSE
		}
	
		Kmin <- Kmin+1
  }

  setEmOptions(emOptions)
   
  for(K in Kmin:Kmax)
  {
   for(xcte in cte)
   {
      ##	Control parameters
      HasExcluded <- FALSE
      HasIncluded <- FALSE
	
	  	K = as.integer(K)
  
      ## Initialise selected set of variables
      S <- rep(TRUE, P)
	
	  	repeat
	  		{
					#####################################
					##	Exclusion step
					#####################################
	
					#	Reference model
					if(!all(S == Smin))
					{
		  			cat("\n>>> In Exclusion with K = ", K, "\n")
	  
		 	 			RefCriteria <- read.or.compute(data2C,
									   K, 
								 S, 
								 xReferenceModel=TRUE, 
								 ReferenceModelsIndex,
								 NberExploredModels,
								 FileNameExploredModels,
								 xcte,
								 header = TRUE)
		
		  			##	Criteria of the reference model
		  			cat("\n>>> K = ", K, "; S = ", c(1:P)[S], "\n")
		  			cat("\n >>> Criteria = ", RefCriteria, "\n")
	  
		  			##	Matrix of criteria of competing models
		  			CompetCriteria <- matrix(NA, ncol=NBER_CRITERIA)
		  			colnames(CompetCriteria) <- NAMES_CRITERIA
	  
		  			for(j in 1:P)
		  			{
							if(S[j] & !Smin[j])
							{
			  				xxS <- S
			  				xxS[j] <- FALSE
			
			  				xxCompetCriteria = read.or.compute(data2C,
											   K, 
									  xxS, 
									  xReferenceModel=FALSE, 
									  ReferenceModelsIndex,
									  NberExploredModels,
									  FileNameExploredModels,
									  xcte,
									  header = TRUE)
			
			  				CompetCriteria = rbind(CompetCriteria, xxCompetCriteria)
			
			  				##	Print xxS
			  				cat("\n>>> K = ", K, "; S[-", j, "] = ", c(1:P)[xxS], "; Criteria = ", xxCompetCriteria, "\n >>> Constant = ", xcte, "\n")
							}
							else CompetCriteria = rbind(CompetCriteria, rep(NA, NBER_CRITERIA))
		  			}
	  
		  			#	Warning: CompetCriteria has been initialized with NAs on the first line
		  			candidateEx = which.min(CompetCriteria[, CRIT]) - 1
		
		  			##	Print exclusion candidate
		  			cat("\n>>> Candidate for Exclusion = ", candidateEx, "\n")
	  
		  			if((ForceExclusion && !HasIncluded) | (RefCriteria[CRIT] >= CompetCriteria[candidateEx + 1, CRIT]))
		  				{
								##	Then exclusion
								S[candidateEx] <- FALSE
								HasExcluded <- TRUE
								RefCriteria <- CompetCriteria[candidateEx+1, ]
		  
								##	Print excluded variable
								cat("\n>>> Variable ", candidateEx, " excluded\n")
		  				}
		  				else HasExcluded = FALSE
					}
					else HasExcluded = FALSE
	  
					###########################################
					##	Inclusion step
					###########################################
	  
					#print("Before Inclusion")
	  
					##	Break if there is no variables to include
					if(!all(S))
					{
		  			#	reference model come from privious Exclusion step
		
		  			##	Debug
		  			cat("\n>>> In Inclusion step with K = ", K, "\n")
		
		  			##	Competing criteria
		  			CompetCriteria <- matrix(NA, ncol=NBER_CRITERIA)
		  			colnames(CompetCriteria) <- NAMES_CRITERIA
	  
		  			for(j in 1:P)
		  			{
							if(!S[j])
							{
			  				xxS <- S
			  				xxS[j] <- TRUE
			  				xxCompetCriteria <- read.or.compute(data2C,
												K, 
									   xxS, 
									   xReferenceModel=FALSE, 
									   ReferenceModelsIndex,
									   NberExploredModels,
									   FileNameExploredModels,
									   xcte,
									   header = TRUE)
			
			  				CompetCriteria = rbind(CompetCriteria, xxCompetCriteria)
			
			  				##	Print xxS
			  				cat("\n>>> K = ", K, "; S[+", j, "] = ", c(1:P)[xxS], "; Criteria = ", xxCompetCriteria, "\n >>> Constant = ", xcte, "\n")
							}
							else CompetCriteria = rbind(CompetCriteria, rep(NA, NBER_CRITERIA))
		  			}
	  
		  			#	Warning: CompetCriteria has been initialized with NAs on the first line
		  			candidateIn <- which.min(CompetCriteria[, CRIT]) - 1
	  
		  			candidateS <- as.logical(S)
		  			candidateS[candidateIn] <- TRUE
		
		  			##	Print inclusion candidate
		  			cat("\n>>> Candidate for Inclusion = ", candidateIn, "\n")
	  
		  			found = isInFile_Rcpp(K, candidateS, FileNameExploredModels, header = TRUE) #list(TrueFalse, line, dim, logLik)
		  			numModel = found[["line"]] + 1
		
		  			if((RefCriteria[CRIT] > CompetCriteria[candidateIn + 1, CRIT]) &  (!ReferenceModelsIndex[numModel]))
		  			{
							##	Then Inclusion
							S[candidateIn] <- TRUE
							HasIncluded <- TRUE
							RefCriteria <- CompetCriteria[candidateIn+1, ]
							ReferenceModelsIndex[numModel] <- TRUE
		  
							##	Print included variable
							cat("\n>>> Variable ", candidateIn, " included\n")
		  			}
		  			else HasIncluded <- FALSE
		
		  			##	Debug
		  			cat("Out Inclusion")
					}
					else HasIncluded <- FALSE
	  
	  
					##	Break repeat
					if(!ForceExclusion)
					{
		  			if(((!HasExcluded) & (!HasIncluded)) | sum(as.integer(S) - as.integer(Smin)) <= 1) break
					}
					if(sum(as.integer(S) - as.integer(Smin)) <= 1) break
				}
			}
  	}
  
    
  	cat("\n\n>>> The explored models are in the file named ", FileNameExploredModels, "\n")
  
  	## Best models indexes
  	##Rcpp
  	ListCriteria = computeCriteriaFromFile_Rcpp(FileNameExploredModels, 1.0, header = TRUE)# From cpp
  
  	bestModelsIndexes = ListCriteria[["bestModelsIndexes"]] + 1
  
  	##Rcpp
  	bestModels = sapply(bestModelsIndexes, function(i) return (readModelAt_Rcpp(FileNameExploredModels, i - 1, header = TRUE)))
  	colnames(bestModels) = getCriteriaNames_Rcpp()
  	rownames(bestModels) = c("N", "P", "K", "S", "logLik", "dim", "entropy")

  	try(initialiseEmSettings())
  
  	invisible(bestModels)
}

## Constant calibration
#'	Data driven calibration of the penalty function.
#'
#'	@title	Data driven calibration of the penalty function
#'
#'	\code{dimJump.R(fileOrData, h = integer(), N = integer(), header = logical())} returns a data driven penalty parameter.
#'	
#'	@param	fileOrData	A character string or a data.frame (see details).
#'	@param	h	An integer defining the sliding window used to find the biggest jump.
#'	@param	N	The size of the sample data.
#'	@param	header	A logical indicating if the file containts header or not.
#'
#'	@return	NULL
#'
#'	@details
#'	This function is a dimension jump version of the slope heuristics for calibration of penalty function.
#'
#'	@seealso \code{\link{model.selection.R}} for final selection.
#'
#'	@keywords Cluster Classification
#'
#'	@author	Wilson Toussile
#'
#'	@references
#'	\itemize{
#'		\item  \href{http://projecteuclid.org/euclid.ejs/1379596773}{Dominique Bontemps and Wilson Toussile (2013)} : Clustering and #'	variable selection for categorical multivariate data. Electronic Journal of Statistics, Volume 7, 2344-2371, ISSN.
#'		\item 	\href{http://www.springerlink.com/content/r7k02q167120xl64/}{Wilson Toussile and Elisabeth Gassiat (2009)} : Variable #'	selection in model-based clustering using multilocus genotype data. Adv Data Anal Classif, Vol 3, number 2, 109-134.
#'	}
#'	
#'	@examples
#'	\dontrun{
#'	data(genotype2_ExploredModels)
#'	outDimJump = dimJump.R(genotype2_ExploredModels, N = 1000, h = 5)
#'	outDimJump[[1]]
#'	}
#'	@export
dimJump.R <- function(fileOrData, h = integer(), N = integer(), header = logical())
{
	## Compatibility of data
	if(!is.character(fileOrData) & (!is.data.frame(fileOrData)))
		stop(gettextf("Argument 'fileOrData' must be a string or a data.frame"))
	
	if(is.character(fileOrData))
	{
		if(!file.exists(fileOrData))
			stop(paste("File ", fileOrData, " does not exists"))
			
		exploredModels = read.table(fileOrData, header = header)
		xlogLik = exploredModels[, "logLik"]
		xdim = exploredModels[, "dim"]
	}
	else
	{
		xlogLik = as.double(fileOrData[, "logLik"])
		xdim = as.integer(fileOrData[, "dim"])
	}
	
	if(N <= 0) stop("The size N of dataset is incorrect")
			
	
	windowSize = 0.075
	maxConstant = log(N)
	minConstant = 0.25
	gridLength = round(log(N)/windowSize) + 1
	vectConstants = as.double(seq(from = windowSize, to = windowSize * gridLength, by = windowSize))
	vectLogLik = double(length(vectConstants))
	vectDim = integer(length(vectConstants))
	
	##Rcpp
	if(!try(selectDimFromData_Rcpp(xlogLik, xdim, vectConstants, vectLogLik, vectDim)))
		stop("Selection of dimension failed")
	
	BeginEnd1 = integer(2)
	BeginEnd2 = integer(2)
	
	##Rcpp
	dimJump_Rcpp(vectDim, h, BeginEnd1, BeginEnd2)
	
	BeginEnd1 = BeginEnd1 + 1
	BeginEnd2 = BeginEnd2 + 1
	
	constant1 = mean(vectConstants[BeginEnd1])
	constant2 = mean(vectConstants[BeginEnd2])
	
	which_max_dim = which.max(xdim)
	a1 = xlogLik[which_max_dim] - constant1*xdim[which_max_dim]
	a2 = xlogLik[which_max_dim] - constant2*xdim[which_max_dim]
	
	dev.new(title = "Penalty parameter", width = 12, height = 6)
	par(mfrow = c(1, 2))
	
	plot(vectConstants, vectDim, pch = 20, col = "grey40", xlab = "Parameter", ylab = "selected dimension")
	abline(v = vectConstants[BeginEnd1[1]], col = "blue", lty = 3)
	abline(v = constant1, col = "blue")
	abline(v = vectConstants[BeginEnd1[2]], col = "blue", lty = 3)
	abline(v = vectConstants[BeginEnd2[1]], col = "red", lty = 3)
	abline(v = constant2, col = "red")
	abline(v = vectConstants[BeginEnd2[2]], col = "red", lty = 3)
	text(constant1, mean(vectDim), paste("Biggest jump = ", constant1), col = "blue")
	text(constant2, 0.50*mean(vectDim), paste("Second biggest jump = ", constant2), col = "red")
	legend("topright", legend = c(paste("Biggest jump at = ", constant1), paste("Second biggest jump at = ", constant2)), lty = 1, col = c("blue", "red"))
	title(main = "Proposed data driven parameter")
	
	plot(xdim, xlogLik, xlab = "dim", ylab = "Log-likelihood", col = "grey40")
	abline(a1, constant1, col = "blue", lwd = 2, lty = 1)
	abline(a2, constant2, col = "red", lwd = 2, lty = 2)
	legend("topleft", legend = c(paste("Slope 1 = ", constant1), paste("Slope 2 = ", constant2)), col = c("blue", "red"), lty = c(1, 2), lwd = 2)
	title(main = "Associated slope")
	
	bounds = cbind(vectConstants[BeginEnd1], vectConstants[BeginEnd2])
	rownames(bounds) = c("lower", "uper")
	colnames(bounds) = c("constant1", "constant2")
	
	invisible(list(constants = c(constant1, constant2), bounds = bounds, dimLogLik = data.frame(dim = vectDim, logLik = vectLogLik)))
}


## Function for model selection
#'	Function for selection of both the number K of clusters and the subset S of clustering variables.
#'	
#'	@title	Function for selection of both the number K of clusters and the subset S of clustering variables.
#'	
#'	\code{model.selection.R(fileOrData, cte, header, lines)} returns a data.frame in which each line defines 
#'	the selected model for a criterion.
#'
#'	@param	fileOrData	A character string or a data.frame (see details). If \code{fileOrData} is a data.frame, 
#'	it must containts at a column named \code{logLik} and another named \code{dim} (see details).
#'	@param	cte	A penalty function parameter. The associated criterion is \code{-log(likelihood)+cte*dim}.
#'	@param	header	Indication of the presence of header in the file.
#'	@param	lines	A vector of integer. If not empty and \code{fileOrData} is the name of a file, only models defined 
#'		in lines in \code{lines} are compared.
#'	
#'	@return	A data.frame of selected models for the proposed penalized criteria.
#'
#'	@author	Wilson Toussile
#'
#'	@references
#'	\itemize{
#'		\item  \href{http://projecteuclid.org/euclid.ejs/1379596773}{Dominique Bontemps and Wilson Toussile (2013)} : Clustering and #'	variable selection for categorical multivariate data. Electronic Journal of Statistics, Volume 7, 2344-2371, ISSN.
#'		\item 	\href{http://www.springerlink.com/content/r7k02q167120xl64/}{Wilson Toussile and Elisabeth Gassiat (2009)} : Variable #'	selection in model-based clustering using multilocus genotype data. Adv Data Anal Classif, Vol 3, number 2, 109-134.
#'	}
#'	
#'	@examples
#'	\dontrun{
#'	data(genotype2_ExploredModels)
#'	outDimJump = dimJump.R(genotype2_ExploredModels, N = 1000, h = 5)
#'	cte1 = outDimJump[[1]][1]
#'	outSlection = model.selection.R(genotype2_ExploredModels, cte = cte1, header = TRUE)
#'	outSlection
#'	}
#'	@export
model.selection.R <- function(fileOrData,
								cte = as.double(1.0), 
								alpha = as.double(2.0),
								header = TRUE, 
								lines = integer())
{
	## Compatibility of data
	if(!is.character(fileOrData) & (!is.data.frame(fileOrData)))
		stop(gettextf("Argument 'fileOrData' must be a character string or a data frame"))
		
	if((alpha < 1.5) | (alpha > 2))
	  stop(gettextf("The slop coefficient 'alpha' must be in [1.5, 2]"))
	  
	cte = alpha * cte
	
	##Rcpp
	NBER_CRITERIA = getNberCriteria_Rcpp()
		
	if(is.character(fileOrData))
	{
		if(!file.exists(fileOrData))
			stop(gettextf("File '%s' does not exists in the working directory", fileOrData))
			
		parList = readModelAt_Rcpp(fileOrData, 0, header)
		P = parList[["P"]]
		#N = parList[["N"]]
			
		lines = as.integer(lines - 1)
			
		vectK = integer(NBER_CRITERIA)
		matS = matrix(as.integer(0), nrow = NBER_CRITERIA, ncol = P)
		vectLogLik = double(NBER_CRITERIA)
		vectDim = integer(NBER_CRITERIA)
		vectEntropy = double(NBER_CRITERIA)
		vectCriteria = double(NBER_CRITERIA)
		vectN = integer(1)
		
		##Rcpp
		try(selectModelFromFile_Rcpp(fileOrData, vectN, vectK, matS, vectLogLik, vectDim, vectEntropy, vectCriteria, cte, header, lines))
		
		outPut = data.frame(N = rep(vectN[1], NBER_CRITERIA), P = rep(P, NBER_CRITERIA), K = vectK, X = matS, logLik = vectLogLik, dim = vectDim, entropy = vectEntropy)
	}
	else
	{
	  ## Revoir
		N = fileOrData[1, 1]
		P = fileOrData[1, 2]
		vectIndexes = integer(NBER_CRITERIA)
		vectCriteria = double(NBER_CRITERIA)
		vectLogLik = fileOrData[, "logLik"]
		vectDim = fileOrData[, "dim"]
		vectEntropy = fileOrData[, "entropy"]
		
		##Rcpp
		try(selectModelFromData_Rcpp(vectLogLik, vectDim, vectEntropy, N, cte, vectIndexes, vectCriteria))
		vectIndexes = vectIndexes + 1
		
		outPut = data.frame(fileOrData[vectIndexes, ], criteria = vectCriteria)
	}
	
	##Rcpp
	rownames(outPut) = getCriteriaNames_Rcpp()
	
	return (outPut)
}


## Documentation of data
#' @title \code{genotype1} is a data frame of genotype data.
#'
#'	@name genotype1
#'
#'	@docType data
#'	
#'	@format
#'	The format is: chr [1:1000, 1:10] "109107" "105101" "106106" ... and the 11-th column containts integers representing 
#'	the prior classification in \code{5} sub-populations.
#'
#'	@description A simulated data frame of genotype data with \code{N = 1000} individuals gegnotyped at \code{P = 10} loci. 
#'	Each string represents two alleles. For example, "109107" represents "109" and "107". The last column 
#'
#'	@author Wilson Toussile
#'
#'	@references
#'	\itemize{
#'		\item  \href{http://projecteuclid.org/euclid.ejs/1379596773}{Dominique Bontemps and Wilson Toussile (2013)} : Clustering and #'	variable selection for categorical multivariate data. Electronic Journal of Statistics, Volume 7, 2344-2371, ISSN.
#'		\item 	\href{http://www.springerlink.com/content/r7k02q167120xl64/}{Wilson Toussile and Elisabeth Gassiat (2009)} : Variable #'	selection in model-based clustering using multilocus genotype data. Adv Data Anal Classif, Vol 3, number 2, 109-134.
#'	}
#'
#'	
#'	@keywords genotype
NULL






