estimates <- function(indiv, csp) {
  # Estimation of dropout values.
  #
  # Parameters: 
  #    indiv: Profile of an individual as a list of locus, each element
  #           containing one or two strings naming the alleles.
  #    csp: Crime-scene profile. Matrix of replicate (rows) vs loci (columns).
  #         Each element is a vector of strings naming the alleles present in
  #         the crime-scene profile. 
  meanrep = array(0, nrow(csp))
  for(rep in 1:nrow(csp))  {
    represented = c()
    for(locus in colnames(csp)) 
      if(!is.null(csp[[rep, locus]])) {
        represented = c(represented, indiv[[locus]] %in% csp[[rep, locus]])
      }
    meanrep[[rep]] <- sum(represented) / length(represented) 
  }
  meanrep <- 1e0 - meanrep / nrow(csp)
  meanrep[meanrep > 0.99] = 0.99
  meanrep[meanrep < 0.01] = 0.01
  meanrep
}

# Best(?) guess for initial arguments. 
initial.arguments <- function(hypothesis, ...) {
  # Best(?) guess for initial arguments. 
  #
  # Parameters: 
  #    hypothesis: Hypothesis for which to guess nuisance paramters.

  hypothesis = add.args.to.hypothesis(hypothesis, ...)
  sanity.check(hypothesis) # makes sure hypothesis has right type.
  # -1 because relative to first.
  nrcont          = max(nrow(hypothesis$dropoutProfs)
                        + hypothesis$nUnknowns - 1, 0)
  locusAdjustment = rep(1, ncol(hypothesis$dropoutProfs))
  dropout         = runif(nrow(hypothesis$cspProfile),min=0.3,max=0.7)
  degradation     = rep( 3e-3, 
                         nrow(hypothesis$dropoutProfs) + hypothesis$nUnknowns )
  rcont           = runif(nrcont, min=0.5, max=1.5)
  dropin          = NULL
  if(hypothesis$doDropin) dropin = 1e-2

  list(locusAdjustment = locusAdjustment,
       power           = -4.35,
       dropout         = dropout, 
       degradation     = degradation,
       rcont           = rcont,
       dropin          = dropin)
}

relistArguments <- function( parameters, hypothesis, fixed=NULL,
                             logDegradation=TRUE, arguments=NULL ) {
  # Remakes arguments from a flat vector into a list.
  if(is.null(arguments)) arguments = initial.arguments(hypothesis) 
  if(!is.null(fixed)) 
       template = arguments[-which(names(arguments) %in% fixed)]
  else template = arguments
  notempty = Filter(function(n) length(n) > 0, template)
  result <- relist(parameters, notempty)
  if(logDegradation && "degradation" %in% names(result))
    result[["degradation"]] = 10^result[["degradation"]]
  if(!is.null(fixed)) result <- append(result, arguments[fixed])
  result
}



upper.bounds = function(arguments, zero=1e-6, logDegradation=FALSE) { 
  # Upper bounds of optimisation function.
  # 
  # Parameters:
  #  arguments: Arguments passed to the optimisation function. Used as a
  #             template.
  #   zero: Some bounds should be given as >, rather than >=. This arguments is
  #         an epsilon to simulate the first case.
  locusAdjustment = rep(1.5, length(arguments$locusAdjustment))
  dropout     = rep(1-zero, length(arguments$dropout))
  degradation = if(logDegradation) { 0-zero } else { 1-zero }
  degradation = rep(degradation, length(arguments$degradation))
  rcont       = rep(100, length(arguments$rcont))
  dropin      = NULL
  if(!is.null(arguments[["dropin"]])) dropin = 10 - zero

  list(locusAdjustment = locusAdjustment,
       power           = -2, 
       dropout         = dropout,
       degradation     = degradation,
       rcont           = rcont,
       dropin          = dropin)[names(arguments)]
}
lower.bounds = function(arguments, zero=1e-6, logDegradation=FALSE) { 
  # Lower bounds of optimisation function.
  # 
  # Parameters:
  #  arguments: Arguments passed to the optimisation function. Used as a
  #             template.
  #   zero: Some bounds should be given as >, rather than >=. This arguments is
  #         an epsilon to simulate the first case.
  #   logDegradation: Wether degradation parameters are entered as exponents of
  #                   10.
  locusAdjustment = rep(0.5, length(arguments$locusAdjustment))
  degradation = if(logDegradation) { -20 } else { 0 }
  degradation = rep(degradation, length(arguments$degradation))
  dropout     = rep(zero, length(arguments$dropout))
  rcont       = rep(zero, length(arguments$rcont))
  dropin      = NULL
  if(!is.null(arguments[["dropin"]])) dropin = zero

  list(locusAdjustment = locusAdjustment,
       power           = -6,
       degradation     = degradation,
       dropout         = dropout,
       rcont           = rcont, 
       dropin          = dropin)[names(arguments)]
}


optimisation.params <- function(hypothesis, verbose=FALSE, fixed=NULL,
                                logObjective=TRUE, logDegradation=TRUE,
                                arguments=NULL, zero=0, throwError=FALSE,
                                withPenalties=TRUE, doLinkage=TRUE, objective=NULL, iterMax=75, likeMatrix=FALSE,...) {
  # Creates the optimisation parameters for optim.
  #
  # optim is the optimisation function from R's stat package.
  # 
  # Parameters:
  #    hypothesis: Hypothesis for which to create the objective function.
  #    verbose: Wether to print each time the likelihood is computed.
  #    fixed: List of arguments to keep fixed.
  #    logObjective: Whether to optimize the log10 of the likelihood.
  #    logDegradation: Whether to input the degradation as 10^d or not.
  #    arguments: Starting arguments. If NULL, then uses initial.arguments to
  #               compute them.
  #    zero: An epsilonic number used to indicate lower and upper bounds which
  #          should be excluded.
  #    throwError: Throw an error if result is infinite

  hypothesis = add.args.to.hypothesis(hypothesis, ...)
  sanity.check(hypothesis) # makes sure hypothesis has right type.
  # If the objective function has not been handed to optimizatio.params,
  # make the objective function
  if(is.null(objective)) objective = create.likelihood.vectors(hypothesis, likeMatrix=likeMatrix,...)

  # Get maximum allele fraction
  maxAF <- getMaxAF(hypothesis) 

  args = arguments
  if(is.null(args)) args = initial.arguments(hypothesis)
  if(logDegradation && "degradation" %in% names(args)) 
    args$degradation = log10(args$degradation)
  
  # Make sure we don't include empty stuff (like rcont sometimes)
  template = Filter(function(n) length(n) > 0, args)

  if(!is.null(fixed)) {
    fixedstuff = args[fixed]
    template = args[-which(names(args) %in% fixed)]
  } else  fixedstuff = NULL

    # Linkage adjustment if brothers
    linkBool = doLinkage&!hypothesis$relationship%in%c(0,1)&hypothesis$hypothesis=="prosecution"
  if(linkBool)
    {
    linkFactor = linkage(hypothesis)
    if(logObjective) linkFactor = log10(linkFactor)
    }

    # result function
  result.function <- function(x) {
    # If a flat vector, reconstruct list. 
    if(typeof(x) == "double")
      x = relistArguments(x, hypothesis, fixed=fixed, arguments=arguments, logDegradation=logDegradation)
    # Otherwise, checks some options.
    else { 
      # Make sure it contains fixed terms.
      if(length(fixedstuff)) x = append(x, fixedstuff)
      # Converts degradation from log10 format.
      if(logDegradation && "degradation" %in% names(x))
        x$degradation = 10^x$degradation
    }

    # If would return negative likelihood skip
    if(hypothesis$doDropin & checkDropin(x, maxAF, hypothesis$nUnknowns+nrow(hypothesis$dropoutProfs)))
	{
	if(logObjective) result = log10(0) else result = 0
	if(verbose) print(result)
	return(-result)
	}

    # Compute objective function.
    result <- do.call(objective, x)
if(likeMatrix==TRUE) return(result)
    # Compute as log if requested, otherwise as product.
    if(withPenalties) {
      if(logObjective)
        result <- sum(log10(result$objectives) + log10(result$penalties))
      else result <- prod(result$objectives) * prod(result$penalties)
    } else if(logObjective) {
      result <- sum(log10(result$objectives))
    } else result <- prod(result$objectives) 
    # Print out if requested.
    if(verbose) {
      # print(unlist(append(x, list(result=result))))
      print(result)
    }
    # If result is infinite, do throw
    if(throwError == TRUE && is.infinite(result)) {
      cat("Objective function is over/underflow: ", result, "\n")
      print(x)
      stop("Objective function over/underflow")
    }
    # If result is infinite make sure it is -Inf
    if(is.infinite(result)|is.na(result)) result = -Inf 
    
    # Linkage adjustment
    if(linkBool)
        {
        if(logObjective==TRUE) 
            {
            result = result+linkFactor
            } else {
            result = result*linkFactor
            }
        }
    # return result
    -result
  }


  
  lower = lower.bounds(args, zero, logDegradation)
  upper = upper.bounds(args, zero, logDegradation)
  lower = lower[names(template)] 
  upper = upper[names(template)] 

 # population size for optimisation
  searchPop = 4*length(unlist(upper))
  # increases for relatedness
  searchPop = round(searchPop * relFactor(hypothesis$relatedness))


  list(#par     = unlist(template), 
       fn      = result.function, 
       lower   = unlist(lower), 
       upper   = unlist(upper),
       #control = list(fnscale=-1, factr=1e7, maxit=500), 
       control = list(strategy=3,NP=searchPop,itermax=iterMax) 
       #method  = "L-BFGS-B",
       #hessian = FALSE )
       )
}


# factor by which to increase the population size for optimisation when taking into account relatedness
# optimisation seems more difficult when related, so need a more thorough search
relFactor = function(relatedness,base=1,fac1=9,fac2=4,fac3=2)
	{
	base + fac1*relatedness[1] + fac2*relatedness[2] + fac3*prod(relatedness)
	}
	

getMaxAF <- function(hypothesis)
	{
	# Get the maximum allele frequency of any alleles in database 
	#
	# Parameters:
	# hypothesis: Hypothesis from which to get the maximum allele fraction
	maxAF <- function (db)
		{
		max(db[,1])
		}
	out <- max(sapply(X=hypothesis$alleleDb, FUN=maxAF))
	return(out)
	}


checkDropin <- function(params, maxAF, nContrib)
	{
	# Check that dropin value will not create negative likelihoods
	#
	# Parameters:
	# 	params: parameters to be checked
	# 	maxAF: maximum allele fraction
	#	nContrib: number of unknowns plus number of contributors subject to dropout
	dropin = params$dropin
	dropout = params$dropout
	# dropin rate is just dropin probability if no reference individual
	if(nContrib==0)
		{
		check = maxAF * dropin
		} else {
		check = maxAF * (dropin*(1-dropout))
		}
	if(any(check>1)) return(TRUE) else return(FALSE)
	}

DEoptimLoop <- function(PARAMS, tolerance=1e-6){
	# Optimises over parameters, while checking for convergence every 50 iterations
	# 	PARAMS: parameters for DEoptim generated by optimisation.params

	# Bogus result to check against at first
	oldresult = 999
	globalBestVal = 999
	globalBestMem = NULL
	bestmemitOut <- NULL
	bestvalitOut <- NULL
	iterOut <- NULL
	nfevalOut <- NULL
	# Set check flag to FALSE
	flag = FALSE
	while(!flag){
		# Run DEoptim for 100 iterations
		results <- do.call(DEoptim,PARAMS)
		# Get results for output
		bestmemitOut <- rbind(bestmemitOut, results$member$bestmemit)
		bestvalitOut <- rbind(bestvalitOut, results$member$bestvalit)
		iterOut <- c(iterOut, results$optim$iter)
		nfevalOut <- c(results$optim$nfeval) 
		# Check if result is the same as the pervious one (50 iterations before)
		tocheck = ifelse(globalBestVal==999,globalBestVal,min(globalBestVal, oldresult))
		condition = abs((results$optim$bestval-tocheck)/tocheck)<tolerance
		if(condition) flag=TRUE
		# If not set the initial population to the last population from the last result
		PARAMS$control$initialpop = results$member$pop
		# set global best result
		if(results$optim$bestval<globalBestVal) 
			{
			globalBestVal = results$optim$bestval
			globalBestMem = results$optim$bestmem	
			}
		# save current result
		oldresult = results$optim$bestval
		}
	# make globalBestVal the output result if it is better than oldresult
	if(globalBestVal<results$optim$bestval)
		{
		results$optim$bestval = globalBestVal
		results$optim$bestmem = globalBestMem
		#print("*** DEoptimLoop suboptimal convergence ***")
		}
	# Update some results to include all steps of loop
	results$member$bestmemit = bestmemitOut
	results$member$bestvalit = bestvalitOut
	results$optim$iter = sum(iterOut)
	results$optim$nfeval = sum(nfevalOut)
	return(results)
	}

get.likely.genotypes = function(hypothesis,params,results,posterior=FALSE,joint=FALSE,prob=ifelse(joint==FALSE,0.1,0.05))
	{
	# Function to return likely genotypes for each individual
	# Finds the marginal probabilities of genotypes, and then
	# returns the most likely genotypes, with their probabilities
	#
	# Parameters:
	#	hypothesis: hypothesis from defence.hypothesis(args)
	#	params: parameters object returned from optimisation.params(hypothesis)
	#	results: results object returned from DEoptimLoop(params)
	#	prob: genotypes with a probability greater than this will be returned

	# transform hypothesis to locus centric
	locusCentricHyp = transform.to.locus.centric(hypothesis)
	# function to return genotype combinations for each locus
	genCombs = function(locusHyp,alleleDb)
		{
		genotypeIndex = likelihood.constructs.per.locus(locusHyp)$genotypes
		genotypes = matrix(rownames(alleleDb)[genotypeIndex],ncol=ncol(genotypeIndex))
		return(genotypes)
		}
	genotypes = mapply(FUN=genCombs,locusHyp=locusCentricHyp,alleleDb=hypothesis$alleleDb)
	# make a new output function to output genotype likelihoods
	newParams = optimisation.params(hypothesis,likeMatrix=TRUE)
	# get genotype likelihoods with parameters returned by previous optimisation (results)
	likes = newParams$fn(results$optim$bestmem)$objectives
	# convert to probabilities
	likes = lapply(likes,FUN=function(x) x/sum(x))
	# if want posterior, return
	if(posterior==TRUE) return(list(genotypes=genotypes,probabilities=likes))
	# if we only want the joint distributions
	if(joint==TRUE) 
		{
		# joint genotypes
		outJoint = mapply(FUN=subGens,genotypes,likes,rotate=TRUE,prob=prob,SIMPLIFY=FALSE)
		# top genotype combination
		topGenotypes = mapply(FUN=function(a,b) a[,which.max(b)],genotypes,likes,SIMPLIFY=TRUE)
		topGenotypes = t(topGenotypes)
		# get top probability
		topProbability = prod(sapply(likes,FUN=max))
		return(list(joint=outJoint,topGenotypes=list(genotype=topGenotypes,probability=topProbability)))
		}
	# function to get marginal probabilities and subset to those with prob>prob
	marginalProbs = function(gens,probs,indiv=1,prob=0.1,top=FALSE)
		{
		marginals = marginal(gens,probs,indiv)
		if(top==TRUE)	
			{
			topMarginals = marginals$genotypes[which.max(marginals$probabilities),,drop=FALSE]
			topProbability = max(marginals$probabilities)		
			return(list(genotype=topMarginals,probability=topProbability))
			}
		subMarginals = subGens(marginals$genotypes,marginals$probabilities,prob)
		return(subMarginals)
		}
	# number of contributors
	ncont = nrow(genotypes[[1]])/2
	# get marginal and subset at every locus for every individual
	out = sapply(1:ncont,FUN=function(x) mapply(FUN=marginalProbs,gens=genotypes,probs=likes,indiv=x,prob=prob,SIMPLIFY=FALSE),simplify=FALSE)
	# order by dropout rate
	rcont = vector(length=ncont)
	rcont[hypothesis$refIndiv] = 1
	rcont[-hypothesis$refIndiv] = results$optim$bestmem[grep("rcont",names(results$optim$bestmem))]
	index = (1:ncont)[rev(order(rcont))]
	out = out[index]
	# get top genotypes for marginals
	topGenotypes = sapply(1:ncont,FUN=function(x) mapply(FUN=marginalProbs,gens=genotypes,probs=likes,indiv=x,prob=prob,top=TRUE,SIMPLIFY=FALSE),simplify=FALSE)
	# get top probabilities for marginals	
	topProbabilities = sapply(1:ncont,FUN=function(y) prod(sapply(topGenotypes[[y]],FUN=function(x) x$probability)),simplify=FALSE)
	topGenotypes = sapply(1:ncont,FUN=function(y) sapply(topGenotypes[[y]],FUN=function(x) x$genotype),simplify=FALSE)
	topGenotypes = topGenotypes[index]	
	topGenotypes = lapply(topGenotypes,FUN=t)
	topProbabilities = topProbabilities[index]
	#return(list(genotypes=genotypes,probabilities=likes))
	return(list(out,topGenotypes=list(genotypes=topGenotypes,probabilities=topProbabilities)))
	}


# function to get marginal probabilities for a single contributor
# for internal use within getLikes()
marginal = function(genotypes,probabilities,indiv=1)
	{
	# genotypes for just 1 contributor
	allGen1 = genotypes[(((2*indiv)-1):(2*indiv)),,drop=FALSE]
	# remove duplicate genotypes
	outGens = unique(t(allGen1))
	# sum probabilities for each identical single genotype
	outProbs = apply(outGens,MARGIN=1,FUN=function(x) sum(probabilities[getMatching(singleGens=allGen1,matchGen=x)]))
	return(list(genotypes=outGens,probabilities=outProbs))
	}

# function to return an index of which genotypes in group match a given genotype
# for internal use in marginal()
getMatching = function(singleGens,matchGen)
	{
	index1 = colSums(apply(singleGens,MARGIN=2,FUN=function(x) x%in%matchGen))
	index2 = colSums(apply(singleGens,MARGIN=2,FUN=function(x) matchGen%in%x))
	outIndex = which((index1+index2)==4)
	return(outIndex)
	}

# function to subset genotypes to those with probability greater than prob
# for internal use within getLikes()
subGens = function(genotypes,probabilities,rotate=FALSE,prob=0.1)
	{
	if(rotate==TRUE) genotypes = t(genotypes)
	genotypes = genotypes[rev(order(probabilities)),,drop=FALSE]
	probabilities = sort(probabilities,decreasing=TRUE)
	index = which(probabilities>prob)
	genotypes = genotypes[index,,drop=FALSE]
	probabilities = probabilities[index]
	return(list(genotypes=genotypes,probabilities=probabilities))
	}

#----------------------------------------------------------------------------------------------------
# In preparation for a GUI built in tcltk, we first want regular estimates and outputs of WoE. Eventually will use tkProgressBar.
# This requires a new function called evaluate() which calls DEoptimLoop(), regularly alternating between Pros and Def so they progress together.
# These values can first be output to the screen, to ensure it all works properly, before integrating into the GUI

# Other changes required will be:
# - remove verbose from optimisation.params()
# Strategy:
# currently DEoptimLoop deals with one hyp, looping until convergence is smaller than a tolerance threshold of 1e-06
# Solution: chop up the tolerance threshold into n tolerance steps between 10 and 1e-5, using a geometric series
# this gives the additional benefit of allowing other DEoptim parmameters to be adjusted at each step. 
# Specifically, we adjust CR, starting at 0.1, ending at 0.7 (again a geometric series), which allows a broader search at the beginning, and a more detailed local serach at the end, a bit like a cooling schedule in a simulated annealing algorithm, wrapped around DEoptim.

geometric.series <- function(start,end,n){ 
	# helper function used by evaluate()
	# start: end: Doesn't matter which way around. Values will always be closer to the smaller one, and thin out towards the larger.
	span <- start/end
	inc <- exp(log(span)/(n-1))
	steps <- c(start,start/cumprod(rep(inc,(n-1))))
return(steps)}


evaluate <- function(P.pars, D.pars, tolerance=1e-5, n.steps=NULL, progBar = TRUE, interim=TRUE, CR.start=0.1, CR.end=0.7, seed.input=NULL){

	# P.pars D.pars: parameter object created by optimisation.params()
	# the smallest convergence threshold (ie for the last step)
	# number of convergence thresholds
	# for each step, run a DEoptimLoop both for P and D, until each converges at that steps accuracy

	begin = Sys.time()

	# set seed
	if(is.null(seed.input)) 
	    {
	    seed.used =  as.integer(Sys.getpid()+as.integer(Sys.time()))
	    } else {
        seed.used = as.integer(seed.input)
	    }
    set.seed(seed.used)

	# combine the outputs outside the loop
	P.bestmemit <- D.bestmemit <- NULL
	P.bestvalit <- D.bestvalit <- NULL
	P.iter <- D.iter <- NULL
	P.nfeval <- D.nfeval <- NULL

	# run first step
	n = 1

	# change DEoptim parameters
	P.pars$control$CR <- CR.start
	D.pars$control$CR <- CR.end
		
	# run DEoptimLoop until convergence at the required step
	P.step <- DEoptimLoop(P.pars,10)
	D.step <- DEoptimLoop(D.pars,10)

	# set global results
	GlobalPval = P.step$optim$bestval
	GlobalDval = D.step$optim$bestval
	GlobalPmem = P.step$optim$bestmem
	GlobalDmem = D.step$optim$bestmem

	# put the Pros outputs into the combined output
	P.bestmemit <- rbind(P.bestmemit, P.step$member$bestmemit)
	P.bestvalit <- rbind(P.bestvalit, P.step$member$bestvalit)
	P.iter <- c(P.iter, P.step$optim$iter)
	P.nfeval <- c(P.nfeval, P.step$optim$nfeval)

	# put the Def outputs into the combined output	
	D.bestmemit <- rbind(D.bestmemit, D.step$member$bestmemit)
	D.bestvalit <- rbind(D.bestvalit, D.step$member$bestvalit)
	D.iter <- c(D.iter, D.step$optim$iter)
	D.nfeval <- c(D.nfeval, D.step$optim$nfeval)	

	# recycle the current pop into the next loop
	P.pars$control$initialpop <- P.step$member$pop
	D.pars$control$initialpop <- D.step$member$pop

	# calculate the weight of evidence. Note, results are + rather than -, so D-P
	WoEtmp <- D.step$optim$bestval - P.step$optim$bestval

	# get standard mean standard deviation of initial optimisation phase
	sdStep = mean(c(sd(P.step$member$bestvalit[1:75]),sd(D.step$member$bestvalit[1:75])))
	# sometimes sd is very low (below 1 e.g. 3locus test)
	# if so set sd to >1 so log2(sd) is positive
	if(sdStep<1) sdStep = 1.5
	# decide how many steps to run
	if(is.null(n.steps)) n.steps = ceiling(log2(sdStep))*4+length(grep("cont",names(D.pars$upper)))

	# retain all the likelihood ratios
	WoE <- numeric(n.steps)
	WoE[n] <- WoEtmp

	# intialise progress bar
	if(progBar) 
		{
		pb = tkProgressBar(title="% Progress",min=0,max=n.steps)
		setTkProgressBar(pb,n,label=paste0("WoE = ",round(WoE[n],2)," bans"))
		} else {
		# percentage progress
		progress <- round(100*n/n.steps)
		# print progress
		print(paste('Weight of evidence:',round(WoE[n],2),'Progress:',progress,'%'))
		}

	# if more than one step, start loop
	if(n.steps>1)
		{
		# adjust tolerance gradually
		tol.steps <- geometric.series(10,tolerance,n.steps)
	
		# adjust DEoptim parameters gradually, so search space is confined more at the end
		CR.steps <- geometric.series(CR.start,CR.end,n.steps)

		for(n in 2:n.steps){

			# generate outputs if interim = TRUE
			if(interim==TRUE){
				interim(P.step,D.step,n,n.steps)
				save(list=ls(),file='interim.RData')
				}
	
			# change DEoptim parameters
			P.pars$control$CR <- CR.steps[n]
			D.pars$control$CR <- CR.steps[n]
		
			# run DEoptimLoop until convergence at the required step
			P.step <- DEoptimLoop(P.pars,tol.steps[n])
			D.step <- DEoptimLoop(D.pars,tol.steps[n])

			# set global results
			if(P.step$optim$bestval<GlobalPval)
				{
				GlobalPval = P.step$optim$bestval
				GlobalPmem = P.step$optim$bestmem
				}
			if(D.step$optim$bestval<GlobalDval)
				{
				GlobalDval = D.step$optim$bestval
				GlobalDmem = D.step$optim$bestmem
				}


			# put the Pros outputs into the combined output
			P.bestmemit <- rbind(P.bestmemit, P.step$member$bestmemit)
			P.bestvalit <- rbind(P.bestvalit, P.step$member$bestvalit)
			P.iter <- c(P.iter, P.step$optim$iter)
			P.nfeval <- c(P.nfeval, P.step$optim$nfeval)

			# put the Def outputs into the combined output	
			D.bestmemit <- rbind(D.bestmemit, D.step$member$bestmemit)
			D.bestvalit <- rbind(D.bestvalit, D.step$member$bestvalit)
			D.iter <- c(D.iter, D.step$optim$iter)
			D.nfeval <- c(D.nfeval, D.step$optim$nfeval)	

			# recycle the current pop into the next loop
			P.pars$control$initialpop <- P.step$member$pop
			D.pars$control$initialpop <- D.step$member$pop

			# calculate the weight of evidence. Note, results are + rather than -, so D-P
			WoE[n] <- D.step$optim$bestval - P.step$optim$bestval

			# progress bar
			if(progBar){
				# update progress bar
				setTkProgressBar(pb,n,label=paste0("WoE = ",round(WoE[n],2)," bans"))
				} 
			# percentage progress
			progress <- round(100*n/n.steps)
			print(paste('Weight of evidence:',round(WoE[n],2),'Progress:',progress,'%'))
			}
		}

	changeFlag=FALSE
	# if global result is better than result from last chunk
	if(P.step$optim$bestval>GlobalPval)
		{
		P.step$optim$bestval = GlobalPval
		P.step$optim$bestmem = GlobalPmem
		print("*** Final prosecution result was not the global optimum - consider re-running optimisation ***")
		changeFlag=TRUE
		}
	if(D.step$optim$bestval>GlobalDval)
		{
		D.step$optim$bestval = GlobalDval
		D.step$optim$bestmem = GlobalDmem
		print("*** Final defence result was not the global optimum - consider re-running optimisation ***")
		changeFlag=TRUE
		}

	if(changeFlag) WoE <- c(WoE,D.step$optim$bestval - P.step$optim$bestval)

	# close progress bar
	if(progBar) close(pb)

	# update final Pros results
	P.results <- P.step
	P.results$member$bestmemit <- P.bestmemit
	P.results$member$bestvalit <- P.bestvalit
	P.results$optim$iter <- P.iter
	P.results$optim$nfeval <- P.nfeval

	# update final Pros results
	D.results <- D.step
	D.results$member$bestmemit <- D.bestmemit
	D.results$member$bestvalit <- D.bestvalit
	D.results$optim$iter <- D.iter
	D.results$optim$nfeval <- D.nfeval

	end = Sys.time()
	runtime = list(elapsed=difftime(end,begin),
			start=begin,
			end=end)

# return all results
return(list(Pros =P.results,Def =D.results, WoE =WoE, seed.used=seed.used,seed.input=seed.input,runtime=runtime))}

# function to output interim report
interim = function(resultsP,resultsD,step,n.steps)
	{
	# title
	suppressWarnings(write.table(c("Interim results"),"Interim.csv",append=FALSE,sep=",",row.names=FALSE,col.names=FALSE))
	# step
	if(step==1)
	    {
	    suppressWarnings(write.table(t(c("Step",step)),"Interim.csv",append=TRUE,sep=",",row.names=FALSE,col.names=FALSE))
	    } else {
        suppressWarnings(write.table(t(c("Step",paste(step,n.steps,sep="/"))),"Interim.csv",append=TRUE,sep=",",row.names=FALSE,col.names=FALSE))
	    }
	# WoE
	suppressWarnings(write.table(t(c("WoE",resultsD$optim$bestval-resultsP$optim$bestval)),"Interim.csv",append=TRUE,sep=",",row.names=FALSE,col.names=FALSE))
	# Lp
	suppressWarnings(write.table(t(c("Lp",-resultsP$optim$bestval)),"Interim.csv",append=TRUE,sep=",",row.names=FALSE,col.names=FALSE))
	# Ld
	suppressWarnings(write.table(t(c("Ld",-resultsD$optim$bestval)),"Interim.csv",append=TRUE,sep=",",row.names=FALSE,col.names=FALSE))
	# paramsP
	suppressWarnings(write.table(t(c("paramsP",resultsP$optim$bestmem)),"Interim.csv",append=TRUE,sep=",",row.names=FALSE))
	# paramsD
	suppressWarnings(write.table(t(c("paramsD",resultsD$optim$bestmem)),"Interim.csv",append=TRUE,sep=",",row.names=FALSE))
	}


evaluate.from.interim <- function(file){

	# check if interim.RData exists
	if(!grep(pattern = 'interim.RData', x=file))stop("File must be the 'interim.RData' file produced by evaluate() when interim=T")
	# define variables to appease the CRAN gods
	n.steps = NULL
	CR.steps = NULL
	tol.steps = NULL
	progBar = NULL
	pb = NULL
	seed.used = NULL
	seed.input=NULL
	# load the 'interim.RData' file produced by evaluate(..., interim=T)
	load(file)

		for(n in n:n.steps){

			# generate outputs if interim = TRUE
			if(interim==TRUE){
				interim(P.step,D.step,n,n.steps)
				save(list=ls(),file='interim.RData')
				}

			# change DEoptim parameters
			P.pars$control$CR <- CR.steps[n]
			D.pars$control$CR <- CR.steps[n]
		
			# run DEoptimLoop until convergence at the required step
			P.step <- DEoptimLoop(P.pars,tol.steps[n])
			D.step <- DEoptimLoop(D.pars,tol.steps[n])

			# set global results
			if(P.step$optim$bestval<GlobalPval)
				{
				GlobalPval = P.step$optim$bestval
				GlobalPmem = P.step$optim$bestmem
				}
			if(D.step$optim$bestval<GlobalDval)
				{
				GlobalDval = D.step$optim$bestval
				GlobalDmem = D.step$optim$bestmem
				}

			# put the Pros outputs into the combined output
			P.bestmemit <- rbind(P.bestmemit, P.step$member$bestmemit)
			P.bestvalit <- rbind(P.bestvalit, P.step$member$bestvalit)
			P.iter <- c(P.iter, P.step$optim$iter)
			P.nfeval <- c(P.nfeval, P.step$optim$nfeval)

			# put the Def outputs into the combined output	
			D.bestmemit <- rbind(D.bestmemit, D.step$member$bestmemit)
			D.bestvalit <- rbind(D.bestvalit, D.step$member$bestvalit)
			D.iter <- c(D.iter, D.step$optim$iter)
			D.nfeval <- c(D.nfeval, D.step$optim$nfeval)	

			# recycle the current pop into the next loop
			P.pars$control$initialpop <- P.step$member$pop
			D.pars$control$initialpop <- D.step$member$pop

			# calculate the weight of evidence. Note, results are + rather than -, so D-P
			WoE[n] <- D.step$optim$bestval - P.step$optim$bestval

			# progress bar
			if(progBar){
				# update progress bar
				setTkProgressBar(pb,n,label=paste0("WoE = ",round(WoE[n],2)," bans"))
				} 
			# percentage progress
			progress <- round(100*n/n.steps)
			print(paste('Weight of evidence:',round(WoE[n],2),'Progress:',progress,'%'))
			}


	changeFlag=FALSE
	# if global result is better than result from last chunk
	if(P.step$optim$bestval>GlobalPval)
		{
		P.step$optim$bestval = GlobalPval
		P.step$optim$bestmem = GlobalPmem
		print("*** Final prosecution result was not the global optimum - consider re-running optimisation ***")
		changeFlag=TRUE
		}
	if(D.step$optim$bestval>GlobalDval)
		{
		D.step$optim$bestval = GlobalDval
		D.step$optim$bestmem = GlobalDmem
		print("*** Final defence result was not the global optimum - consider re-running optimisation ***")
		changeFlag=TRUE
		}

	if(changeFlag) WoE <- c(WoE,D.step$optim$bestval - P.step$optim$bestval)

	# close progress bar
	if(progBar) close(pb)

	# update final Pros results
	P.results <- P.step
	P.results$member$bestmemit <- P.bestmemit
	P.results$member$bestvalit <- P.bestvalit
	P.results$optim$iter <- P.iter
	P.results$optim$nfeval <- P.nfeval

	# update final Pros results
	D.results <- D.step
	D.results$member$bestmemit <- D.bestmemit
	D.results$member$bestvalit <- D.bestvalit
	D.results$optim$iter <- D.iter
	D.results$optim$nfeval <- D.nfeval

# return all results
return(list(Pros =P.results,Def =D.results, WoE =WoE, seed.used=seed.used,seed.input=seed.input))}


