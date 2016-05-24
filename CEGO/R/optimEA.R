
###################################################################################
#' Evolutionary Algorithm for Combinatorial Optimization
#' 
#' A basic implementation of a simple Evolutionary Algorithm for Combinatorial Optimization. Default evolutionary operators
#' aim at permutation optimization problems.
#'
#' @param x Optional start individual(s) as a list. If NULL (default), \code{creationFunction} (in \code{control} list) is used to create initial design. 
#' If \code{x} has less individuals than the population size, creationFunction will fill up the rest.
#' @param fun target function to be minimized
#' @param control (list), with the options\cr
#' \code{budget} The limit on number of target function evaluations (stopping criterion) (default: 1000)\cr
#' \code{popsize} Population size (default: 100)\cr
#' \code{generations} Number of generations (stopping criterion) (default: Inf)\cr
#' \code{targetY} Target function value (stopping criterion) (default: -Inf)\cr
#' \code{vectorized} Boolean. Defines whether target function is vectorized (takes a list of solutions as argument) or not (takes single solution as argument). Default: FALSE\cr
#' \code{verbosity} Level of text output during run. Defaults to 0, no output.\cr
#' \code{plotting} Plot optimization progress during run (TRUE) or not (FALSE). Default is FALSE.\cr
#' \code{archive} Whether to keep all candidate solutions and their fitness in an archive (TRUE) or not (FALSE). Default is TRUE.\cr
#' \code{recombinationFunction} Function that performs recombination, default: \code{\link{recombinationPermutationCycleCrossover}}, which is cycle crossover for permutations.\cr
#' \code{recombinationParameters} Parameter list for recombination (e.g., recombinationParameters$recombinationRate => recombination rate, defaults to 0.5). List is passed to recombinationFunction. \cr
#' \code{mutationFunction} Function that performs mutation, default: \code{\link{mutationPermutationSwap}}, which is swap mutation for permutations.\cr
#' \code{mutationParameters} Parameter list for mutation (e.g., mutationParameters$mutationRate => mutation rate). List is passed to mutationFunction. Default: empty list.\cr
#' \code{selection} Selection process: "tournament" (default) or "truncation"\cr
#' \code{tournamentSize} Tournament size (default: 2)\cr
#' \code{tournamentProbability} Tournament probability (default: 0.9)\cr
#' \code{localSearchFunction} If specified, this function is used for a local search step. Default is NULL. \cr
#' \code{localSearchRate} Specifies on what fraction of the population local search is applied. Default is zero. Maximum is 1 (100 percent).\cr
#' \code{localSearchSettings} List of settings passed to the local search function control parameter.\cr
#' \code{stoppingCriterionFunction} Custom additional stopping criterion. Function evaluated on the population, receiving all individuals (list) and their fitness (vector). If the result is FALSE, the algorithm stops.\cr
#' \code{verbosity} >0 for text output.\cr
#' \code{creationFunction} Function to create individuals/solutions in search space. Default is a function that creates random permutations of length 6 \cr
#' \code{duplicateFunction} Function that evaluates a list of solutions for duplicates. Default is the \code{duplicated} function. \cr
#' \code{duplicateRemoval} If set to "all" and archiving is on, new individuals are compared against the complete archive to avoid duplicates. With "population", this is limited to the current population.
#'
#' @return a list:\cr 	
#' \code{xbest} best solution found\cr
#' \code{ybest} fitness of the best solution\cr
#' \code{x} history of all evaluated solutions\cr
#' \code{y} corresponding target function values f(x)\cr
#' \code{count} number of performed target function evaluations 
#' \code{message} Termination message: Which stopping criterion was reached.
#' \code{population} Last population
#' \code{fitness} Fitness of last population
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
#' res <- optimEA(,lF,list(creationFunction=cF,mutationFunction=mF,recombinationFunction=rF,
#'		popsize=15,budget=100,targetY=0,verbosity=1,
#'		vectorized=TRUE)) ##target function is "vectorized", expects list as input
#' res$xbest 
#'
#' @seealso \code{\link{optimCEGO}}, \code{\link{optimRS}}, \code{\link{optim2Opt}} 
#' 
#' @export
###################################################################################
optimEA <- function(x=NULL,fun,control=list()){ 
	#default controls:
	con<-list(budget = 1000 #default controls:
			 , popsize = 100
			 , generations = Inf
			 , targetY = -Inf
			 , vectorized=FALSE
			 , creationFunction = solutionFunctionGeneratorPermutation(6)
			 , recombinationParameters = list(recombinationRate=0.5)
			 , mutationParameters = list()
			 , mutationFunction = mutationPermutationSwap
			 , recombinationFunction = recombinationPermutationCycleCrossover
			 , selection = "tournament" #or "truncation
			 , tournamentSize = 2 
			 , tournamentProbability = 0.9
			 , localSearchFunction = NULL
			 , localSearchRate = 0
			 , localSearchSettings = list()
       , duplicateRemoval = "all"
			 , duplicateFunction = duplicated
       , archive = TRUE
			 , stoppingCriterionFunction = NULL
			 , verbosity = 0 
       , plotting = FALSE			 
			 );
	con$recombinationParameters[names(control$recombinationParameters)] <- control$recombinationParameters
	con$mutationParameters[names(control$mutationParameters)] <- control$mutationParameters
	control$recombinationParameters <- con$recombinationParameters
	control$mutationParameters <- con$mutationParameters
	con[names(control)] <- control;
	control<-con;

  archive <- control$archive
  duplicateRemoval <- control$duplicateRemoval	
	duplicateFunction <- control$duplicateFunction
	creationFunction <- control$creationFunction	
	budget <- control$budget
	vectorized <- control$vectorized
	popsize <- control$popsize
	generations <- control$generations
	targetY <- control$targetY
	tournamentSize <- control$tournamentSize
	recombinationParameters <- control$recombinationParameters
	recombinationRate <- recombinationParameters$recombinationRate
	recombinationFunction <- control$recombinationFunction
	mutationParameters <-  control$mutationParameters
	mutationFunction <- control$mutationFunction
	localSearchFunction <- control$localSearchFunction
	localSearchRate <- control$localSearchRate
	localSearchSettings <- control$localSearchSettings
	localSearchSettings$vectorized <- vectorized
	selection <- control$selection
	stoppingCriterionFunction <- control$stoppingCriterionFunction
	verbosity <- control$verbosity
	plotting <- control$plotting
	tournamentProbability <-  control$tournamentProbability #probability in the tournament selection
	tournamentSize=max(tournamentSize,1)
	
	## Create initial population
	population <- designRandom(x,creationFunction,popsize)
	
	if(vectorized) 
		fitness <- fun(population)
	else
		fitness <- unlist(lapply(population,fun))
		
	count <- popsize	
	gen <- 1
	
  fitbest <- min(fitness,na.rm=TRUE)
  xbest <- population[[which.min(fitness)]]
  if(archive){
    fithist <- fitness
    xhist <- population
  }
	besthist <- fitbest
	run <- TRUE
	while((count < budget) & (gen < generations) & (fitbest > targetY) & (run)){
		gen <- gen+1
		#recombine
		if(selection == "tournament"){ #tournament selection
			popindex <- tournamentSelection(fitness,tournamentSize,tournamentProbability,max(popsize*recombinationRate,2))
		}else{ #truncation selection
			popindex <- order(fitness)[1:max(popsize*recombinationRate,2)]
		}
		offspring <- recombinationFunction(population[popindex],recombinationParameters)
		#mutate
		offspring <- mutationFunction(offspring,mutationParameters)
		#optional local search
		if(!is.null(localSearchFunction) & localSearchRate>0){
			if(localSearchRate<1){
				subsetsize = ceiling(length(offspring)*localSearchRate)
				offspringsubset <- sample(length(offspring),subsetsize)
			}else{
				offspringsubset = 1:length(offspring)
			}	
			for(i in offspringsubset){
				res <- localSearchFunction(x=offspring[i],fun=fun,control=localSearchSettings) 
				offspring[[i]] <- res$xbest #todo: local search already evaluations xbest, is reevaluated in population.
				count <- count + res$count #add local search counted evaluations to evaluation counter of main loop.
			}
		}
		## remove offspring which violate the budget
		offspring <- offspring[1:min(budget-count,length(offspring))]
		## append offspring to population, but remove duplicates first. duplicates are replaced by random, unique solutions.		
    if(duplicateRemoval=="all" & archive) #do not ever evaluate duplicates (compare with archive of solutions)
      offspring <- removeDuplicatesOffspring(xhist,offspring, creationFunction,duplicateFunction)
    else if(duplicateRemoval=="population") #only care for duplicates in current population
      offspring <- removeDuplicatesOffspring(population,offspring, creationFunction,duplicateFunction)
		population <- c(population,  offspring)
		#if any new were created:	
		if(length(population)>length(fitness)){
			#evaluate
			if(vectorized)
				newfit <- fun(offspring)
			else
				newfit <- unlist(lapply(offspring,fun))
      if(archive){
        xhist <- append(xhist,offspring)  
        fithist <-  c(fithist, newfit)
      }  
			fitness <- c(fitness, newfit) #evaluate the new individuals after recombination and mutation
			#update count
			count=count+ length(fitness)-popsize			
			#tournament selection 
			if(selection == "tournament"){ #tournament selection
				popindex <- tournamentSelection(fitness,tournamentSize,tournamentProbability,popsize) #todo should it really be possible to select the same individual several times for the next generation?
			}else{ # truncation selection
				popindex <- order(fitness)[1:popsize]
			}
			population <- population[popindex]
			fitness <- fitness[popindex]			
		}	
    popbest <- min(fitness,na.rm=TRUE)
    if(popbest < fitbest){
      fitbest <- popbest
      xbest <- population[[which.min(fitness)]]
    }   
		if(!is.null(stoppingCriterionFunction)) # calculation of additional stopping criteria
			run <- stoppingCriterionFunction(population,fitness)
		if(plotting){
      besthist <- c(besthist,fitbest)
			plot(besthist,type="l")
		}
    if(verbosity > 0){
    	print(paste("Generations: ",gen," Evaluations: ",count, "Best Fitness: ", min(fitness,na.rm=TRUE)))
    }
	}
	#stopping criteria information for user:
	msg <- "Termination message:"
	if(!run) #success
		msg=paste(msg,"Custom stopping criterion satisfied.")
	if(min(fithist,na.rm=TRUE) <= targetY) 
		msg=paste(msg,"Successfully achieved target fitness.")		
	else if(count >= budget) #budget exceeded
		msg=paste(msg,"Target function evaluation budget depleted.")		
	else if(gen >= generations) #generation limit exceeded
		msg=paste(msg,"Number of generations limit reached.")		
  
  if(archive)
    return(list(xbest=xbest,ybest=fitbest,x=xhist,y=fithist, count=count, message=msg, population=population, fitness=fitness))
  else
    return(list(xbest=xbest,ybest=fitbest,count=count, message=msg, population=population, fitness=fitness)) 
}
