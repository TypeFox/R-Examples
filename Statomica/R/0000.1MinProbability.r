#.packageName <- "LocalFDR"
# empiricalBayes, by Zahra Montazeri and David R. Bickel.

# .packageName <- "empiricalBayes"

default.alternative.prob.threshold <- function(nfeatures)
{
	if(missing(nfeatures) || nfeatures < 1e5)
		c(seq(0.1,0.8,by=0.1),0.85,0.90,0.95,0.96,0.97,0.98,0.99)
	else
		c(0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.98)
}

default.tolerance<-function(x)
{
	th<-sort(x)
	differs<-abs(th[2:length(th)]-th[1:(length(th)-1)])
	min(differs)/2
}
find.alternative.prob.threshold<-function(p.values, alternative.prob.threshold=default.alternative.prob.threshold(), marginal.probability= NULL, max.iteration=10,tolerance=default.tolerance(default.alternative.prob.threshold()), plot.relative.gain = FALSE, call.browser=FALSE, verbose = TRUE)
{
# Returns a numeric vector of length equal to the length of pvalues, indicating the minimum probability that each alternative hypothesis can have to be considered true (conditional on the p.value).
if(!is.numeric(p.values)) stop('p.values should be numeric')
	  is.probability <- function(prob){is.finite(prob) & prob >= 0 & prob <= 1}
	  pvals <- ifelse(is.na(p.values), 1, p.values)
	if(!all(is.probability(pvals))) stop('Each element of p.values must be between 0 and 1')
	if(!is.null(marginal.probability)) 
  {
  if(length(marginal.probability) != 1 || !is.finite(marginal.probability)) stop('marginal.probability should be a single number')
  if(!is.probability(marginal.probability)) stop('marginal.probability must be between 0 and 1')
  }
  if(!all(is.probability(alternative.prob.threshold))) stop('min.probability must be between 0 and 1')
	pi1 <- if(is.null(marginal.probability))
	{
		marginal.probability(p.values=p.values, min.probability=0.5, max.iteration=max.iteration, tolerance=tolerance, verbose = verbose)
	}
	else
		marginal.probability
	if(verbose)
		cat('Using marginal.probability estimate of', pi1, '\n')
	TFvector <- function(j)
	{
		if(verbose)
			message("Computing alternative.probable for threshold ", j, " of ", length(alternative.prob.threshold), " on ", date(), ".")
		alternative.probable(p.values, min.probability=alternative.prob.threshold[j],  marginal.probability = pi1, max.iteration=max.iteration,tolerance=tolerance, plot.relative.gain = plot.relative.gain, call.browser=call.browser)
	}
	TF<-sapply(1:length(alternative.prob.threshold), TFvector)
	p1<- c(0,length(p.values))
	p1scalar<- function(i)
	{
		sca <- if(is.na(TF[i]))
			TF[i]
		else
		{
			if(sum(TF[i,])==0)
				0
			else
				alternative.prob.threshold[sum(TF[i,])]
		}
		stopifnot(length(sca) == 1)
		sca
	}
	p1 <- sapply(1:length(p.values),p1scalar)
# if(!is.null(names(pvalue)))  
 names(p1)<-names(p.values)
 p1

}

priorFDR <- function(p.values, ...) # for users
{
	1 - marginal.probability(p.values = p.values, ...)
}

localFDR <- function(p.values, threshold, prior.fdr, tolerance, ...) # main function for users; returns vector of local false discovery rate estimates
{
	if(missing(threshold))
	{
		threshold <- 1 - default.alternative.prob.threshold(nfeatures = length(p.values))
		message("threshold = ")
		print(threshold)
	}
	p <- 1 - threshold
	get.prior.fdr.tolerance <- function(threshold)
	{
		th <- sort(threshold)
		diffs <- abs(th[2:length(th)] - th[1:(length(th) - 1)])
		min(diffs) / 2
	}
	if(missing(tolerance))
	{
		tolerance <- get.prior.fdr.tolerance(threshold = threshold)
		message("tolerance = ", tolerance)
	}
	marginal.probability <- if(missing(prior.fdr))
		NULL
	else
		1 - prior.fdr
	1 - find.alternative.prob.threshold(p.values = p.values, alternative.prob.threshold = p, marginal.probability = marginal.probability, tolerance = tolerance, ...)
}

