# Calculate the Gini index of <x>
Gini <- function(x) 
{
    n <- length(x)
    x <- sort(x)
    1/(n-1)*(n+1-2*(sum((n+1-seq_len(n))*x)/sum(x)))
}

# Compare the Hamming distances of the successor states of states and their perturbed copies.
testTransitionRobustness <- function(network, accumulate=TRUE, params=list())
{
  res <- do.call("perturbTrajectories", c(list(network=network, measure="hamming"), params))
  if (accumulate)
    return(res$value)
  else
    return(res$stat)
}

# Test attractor robustness by searching the original attractor
# in <params$numSamples> perturbed copies of <network> or by perturbing the state trajectories
# and checking whether the attractors change.
testAttractorRobustness <- function(network, accumulate=TRUE, params=list())
{
  origAttrs <- getAttractors(network,canonical=TRUE)

  # decode parameters
  if (is.null(params$perturb))
    perturb <- "functions"
  else
    perturb <- params$perturb
    
  
  if (params$perturb == "trajectories")
  {
    params$perturb <- NULL
    
    if (!is.null(params$copies))
    {
      params$numSamples <- params$copies
      params$copies <- NULL
    }
    res <- do.call("perturbTrajectories", c(list(network=network, measure="attractor"), params))
    if (accumulate)
      return(res$value)
    else
      return(res$stat)
  }  
  else 
  {
    if (is.null(params$method))
      method <- "bitflip"
    else
      method <- params$method
      
    if (is.null(params$maxNumBits))
      maxNumBits <- 1
    else
      maxNumBits <- params$maxNumBits
      
    if (is.null(params$numStates))
      numStates <- max(1,2^length(network$genes)/100)
    else
      numStates <- params$numStates
      
    if (is.null(params$simplify))
      simplify <- (perturb[1] == "states")
    else
      simplify <- params$simplify
    
    if (is.null(params$readableFunctions))
      readableFunctions <- FALSE
    else
      readableFunctions <- params$readableFunctions
      
    if (is.null(params$excludeFixed))
      excludeFixed <- TRUE
    else
      excludeFixed <- params$excludeFixed

    if (!is.null(params$numSamples))
      copies <- params$numSamples
    else    
    if (is.null(params$copies))
      copies <- 100
    else
      copies <- params$copies
    
    perturbationResults <- unlist(sapply(seq_len(copies),function(copy)
      {
        # get attractors of perturbed network
        perturbedAttrs <- getAttractors(perturbNetwork(network,
                                                       perturb=perturb,
                                                       method=method,
                                                       maxNumBits=maxNumBits,
                                                       numStates=numStates,
                                                       simplify=simplify,
                                                       readableFunctions=readableFunctions,
                                                       excludeFixed=excludeFixed))
        
        # try to find original attractors in perturbed network                                               
	      attractorIndices <- sapply(origAttrs$attractors,function(attractor1)
				  {
					  index <- which(sapply(perturbedAttrs$attractors,function(attractor2)
						  {
							  identical(attractor1,attractor2)
						  }))
					  if (length(index) == 0)
						  NA
					  else
						  index
				  })
	      return(sum(!is.na(attractorIndices)))
	    }))

	  if (accumulate)
	    # return overall percentage of found attractors
    	return(sum(perturbationResults)/(length(origAttrs$attractors) * copies) * 100)
    else
      # return percentage of found attractors for each run
      return(perturbationResults/(length(origAttrs$attractors)) * 100)
  }
}

# Calculate the in-degrees of states in the network.
# If <accumulate> is true, the in-degrees are accumulated
# using the Gini coefficient
testIndegree <- function(network,accumulate=TRUE,params)
{
  attr <- getAttractors(network)
  graph <- plotStateGraph(attr,plotIt=FALSE)
  
  # calculate in-degree using igraph
  if (accumulate)
    # accumulate using Gini index
    return(Gini(degree(graph,mode="in",loops=TRUE)))
  else
    # return the raw degrees
    return(degree(graph,mode="in",loops=TRUE))
}

# Calculate the Kullback-Leibler distance of
# two distributions <x> and <y>.
# <bins> is the number of bins used for discretization.
# <minVal> is the minimum value to be used instead of zero
kullbackLeiblerDistance <- function(x,y,bins=list(),minVal=0.00001)
{
  x <- cut(x,breaks=bins,include.lowest=T,right=F)
  y <- cut(y,breaks=bins,include.lowest=T,right=F)

  tx <- table(x)
  ty <- table(y)
  
  tx <- tx/length(x)
  ty <- ty/length(y)
  
  tx[tx < minVal] <- minVal
  ty[ty < minVal] <- minVal
    
  return(sum(tx*(log(tx/ty))))  
}

# Generic function to test properties of <network> against random networks.
# <numRandomNets> specifies the number of random networks to generate.
# <testFunction> is the name of a function that returns a distribution of test values or a test statistic for each network,
# which receives <testFunctionParams> as optional parameters.
# If <accumulation> is "characteristic", the test function must return a characteristic/test statistic value for each value, and
# a histogram of these value is plotted with a line for the original network.
# If <accumulation> is "kullback_leibler", a histogram of the Kullback-Leibler distances of the test value distribution for the original network
# and the random networks is plotted.
# <sign.level> is the desired significance level for <network> in comparison to the random networks.
# <functionGeneration>,<simplify>,<noIrrelevantGenes>,<validationFunction>,<failureInterations,
# <d_lattice>,<zeroBias> are the corresponding parameters of generateRandomNKNetwork().
# <title> is the title of the plot, <xlab> is its x axis caption, <breaks> is the corresponding histogram parameter, and ... supplies further
# graphical parameters
testNetworkProperties <- function(network, numRandomNets=100, testFunction="testIndegree",
                                  testFunctionParams=list(),
                                  accumulation=c("characteristic","kullback_leibler"),
                                  alternative=c("greater","less"),
                                  sign.level=0.05, drawSignificanceLevel=TRUE,
                                  klBins,klMinVal=0.00001,
                                  linkage=c("uniform","lattice"),
                                  functionGeneration=c("uniform","biased"),
                                  validationFunction, failureIterations=10000,
                                  simplify=FALSE, noIrrelevantGenes=TRUE,
                                  d_lattice=1, zeroBias=0.5, 
                                  title="", xlab, xlim, breaks=30, 
                                  ...)
{
  stopifnot(inherits(network,"BooleanNetwork") || inherits(network,"SymbolicBooleanNetwork"))

  if (is.character(testFunction))
    testFunctionName <- testFunction
  else
    testFunctionName <- ""
  
  testFunction <- match.fun(testFunction)
  
  accumulate <- (match.arg(accumulation) == "characteristic")
  origResult <- testFunction(network,accumulate,testFunctionParams)

  numGenes <- length(network$interactions)
  
  if (inherits(network,"SymbolicBooleanNetwork"))
    inputGenes <- sapply(network$interactions,function(interaction)length(getInputs(interaction)))
  else
    inputGenes <- sapply(network$interactions,function(interaction)length(interaction$input))
  
  if (missing(validationFunction))
    validationFunction <- NULL
      
  randomResults <- lapply(seq_len(numRandomNets),function(i)
                   {      
                      randomNet <- generateRandomNKNetwork(n=numGenes,
                                                           k=inputGenes,
                                                           topology="fixed",
                                                           linkage=linkage,
                                                           functionGeneration=functionGeneration,
                                                           validationFunction=validationFunction,
                                                           failureIterations=failureIterations,
                                                           simplify=simplify,
                                                           noIrrelevantGenes=noIrrelevantGenes,
                                                           d_lattice=d_lattice,
                                                           zeroBias=zeroBias)
                      randomRes <- testFunction(randomNet,accumulate,testFunctionParams)
                      return(randomRes)
                   })
  if (accumulate)
    randomResults <- unlist(randomResults)                 
  
  args <- list(...)
                   
  res <- switch(match.arg(accumulation,c("characteristic","kullback_leibler")),
    characteristic = {
                  # get one value for each random network, and plot a histogram of these values
                  if (missing(xlab))
                  {                    
                    xlab <- switch(testFunctionName,
                      testIndegree = "Gini index of state in-degrees",
                      testAttractorRobustness = "% of identical attractors",
                      testTransitionRobustness = "Normalized Hamming distance",
                      "accumulated results"
                    )
                  }
                  if (missing(xlim))
                  {                    
                    #xlim <- switch(testFunctionName,
                    #  testIndegree = c(0,1),
                    #  testAttractorRobustness = c(0,100),
                    #  c(min(c(origResult,randomResults)),
                    #        max(c(origResult,randomResults)))
                    #)
                    xlim <- range(c(origResult, randomResults))
                  }
                  
                  alternative <- match.arg(alternative, c("greater","less"))
                  
                  # calculate p-value
                  
                  if (alternative == "greater")
                    pval <- sum(randomResults < origResult) / length(randomResults)
                  else                    
                    pval <- sum(randomResults > origResult) / length(randomResults)
                    
                  # plot histogram
                  if (testFunctionName == "testIndegree" | testFunctionName == "testAttractorRobustness")
                  {
                    r <- hist(randomResults,xlim=xlim,xlab=xlab,main=title,xaxt="n",...)
                    axis(side=1,at=seq(xlim[1],xlim[2],length.out=11))
                  }
                  else
                  {
                    # plot with default axis
                    r <- hist(randomResults,xlim=xlim,xlab=xlab,main=title,...)
                  }
                            
                  # plot result for original network          
                  abline(v=origResult,col="red")
                  
                  if (alternative == "greater")
                    text(x=origResult,pos=2,y=max(r$counts)*0.75,
                         labels=paste("> ",round(pval * 100),"%\nof random results",sep=""), 
                         col="red",cex=0.75)
                  else
                    text(x=origResult,pos=4,y=max(r$counts)*0.75,
                         labels=paste("< ",round(pval * 100),"%\nof random results",sep=""), 
                         col="red",cex=0.75)
                       
                  if (drawSignificanceLevel)
                  # plot line for significance level
                  {
                    if (alternative == "greater")
                    {
                      quant <- quantile(randomResults,1.0-sign.level)
                      abline(v=quant,col="blue")
                      text(x=quant,pos=2,y=max(r$counts)*0.85,
                         labels=paste((1.0-sign.level) * 100,"% quantile",sep=""), 
                         col="blue",cex=0.75)                    
                    }
                    else
                    {
                      quant <- quantile(randomResults,sign.level)
                      abline(v=quant,col="blue")
                      text(x=quant,pos=4,y=max(r$counts)*0.85,
                         labels=paste(sign.level * 100,"% quantile",sep=""), 
                         col="blue",cex=0.75)
                    }
                  }
                  list(hist=r,pval=1.0-pval,significant=(1.0-pval<=sign.level))
                },
    kullback_leibler = 
                {
                  # a distribution of values is returned for each network,
                  # plot the Kullback-Leibler distances to the original network
                  
                  if (missing(xlab))
                    xlab <- "Kullback-Leibler distance"
                    
                  if (missing(klBins))
                  {
                     bins <- unique(c(origResult,unlist(randomResults)))
                     bins <- c(bins,max(bins) + 1)
                  }
                  else
                  {
                    bins <- unique(c(origResult,unlist(randomResults)))

                    if (klBins < length(bins))
                      bins <- seq(min(bins),max(bins),length.out=klBins+1)
                    else
                      bins <- c(bins,max(bins) + 1)
                  }
                    
                  vals <- sapply(randomResults,function(results)
                                  kullbackLeiblerDistance(origResult,results,bins=bins,minVal=klMinVal))

                  r <- hist(vals,xlab=xlab,main=title,breaks=breaks,...)
                  list(hist=r)
                },
    stop("'accumulation' must be one of \"characteristic\",\"kullback_leibler\""))   
  return(res)
}
