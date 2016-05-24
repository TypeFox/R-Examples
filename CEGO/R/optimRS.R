
###################################################################################
#' Combinatorial Random Search
#' 
#' Random Search for mixed or combinatorial optimization. Solutions are generated completely at random.
#'
#' @param x Optional set of solution(s) as a list, which are added to the randomly generated solutions and are also evaluated with the target function.
#' @param fun target function to be minimized
#' @param control (list), with the options\cr
#' \code{budget} The limit on number of target function evaluations (stopping criterion) (default: 100)\cr
#' \code{vectorized} Boolean. Defines whether target function is vectorized (takes a list of solutions as argument) or not (takes single solution as argument). Default: FALSE\cr
#' \code{creationFunction} Function to create individuals/solutions in search space. Default is a function that creates random permutations of length 6
#'
#' @return a list:\cr 	
#' \code{xbest} best solution found\cr
#' \code{ybest} fitness of the best solution\cr
#' \code{x} history of all evaluated solutions\cr
#' \code{y} corresponding target function values f(x)\cr
#' \code{count} number of performed target function evaluations 
#' 
#' @examples
#' seed=0
#' glgseed=1
#' #distance
#' dF <- distancePermutationHamming
#' #mutation
#' mF <- mutationPermutationSwap
#' #recombination
#' rF <-  recombinationPermutationCycleCrossover 
#' #creation
#' cF <- function()sample(5)
#' #objective function
#' lF <- landscapeGeneratorUNI(1:5,dF)
#' #start optimization
#' set.seed(seed)
#' res <- optimRS(,lF,list(creationFunction=cF,budget=100,
#'		vectorized=TRUE)) ##target function is "vectorized", expects list as input
#' res$xbest 
#'
#' @seealso \code{\link{optimCEGO}}, \code{\link{optimEA}} 
#' 
#' @export
###################################################################################
optimRS <- function(x=NULL,fun,control=list()){
	con<-list( budget=100
						,vectorized=FALSE
						,creationFunction = solutionFunctionGeneratorPermutation(6)
			 )
	con[names(control)] <- control
	control<-con
	
  budget <- control$budget
	vectorized <- control$vectorized
 	creationFunction <- control$creationFunction	
	
	## Create random solutions without duplicates, filled up with x
	x <- designRandom(x,creationFunction,budget)
  
  #evaluate
  if(vectorized) 
		y <- fun(x)
	else
		y <- unlist(lapply(x,fun))
  
  #best value found:  
	j <- which.min(y)
  
  #return
	list(xbest=x[[j]],ybest=y[j],x=x,y=y, count=budget)
}

