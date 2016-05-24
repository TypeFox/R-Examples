# HighProbability, by David R. Bickel and Zahra Montazeri.

# .packageName <- "HighProbability"

get.marginal.probability.tolerance <- function()
{
	0.005
}


alternative.beneficial <- function(p.values, cost.to.benefit = 1, marginal.probability = NULL,  max.iteration=10, tolerance=get.marginal.probability.tolerance(), plot.relative.gain = FALSE)
# Returns a boolean vector of length equal to the length of p.values, indicating whether each alternative hypothesis is beneficially considered true, i.e., whether its acceptance is optimal given the cost.to.benefit ratio.
{
  if(!is.finite(cost.to.benefit) || length(cost.to.benefit) != 1 || cost.to.benefit <= 0) stop('cost.to.benefit must be a number greater than 0')
  min.probability <- 1 - 1 / (1 + cost.to.benefit)
  cat(paste('Calling alternative.probable with min.probability =', min.probability), '\n')
  alternative.probable(p.values = p.values, min.probability = min.probability, marginal.probability = marginal.probability,  max.iteration=max.iteration, tolerance=tolerance ,plot.relative.gain = plot.relative.gain, call.browser=FALSE)
}

alternative.probable<- function(p.values, min.probability = 0.5, marginal.probability = NULL,  max.iteration=10, tolerance=get.marginal.probability.tolerance(), plot.relative.gain = FALSE, call.browser=FALSE)
# Returns a boolean vector of length equal to the length of p.values, indicating whether each alternative hypothesis is probable, i.e., whether it has a probability (conditional on the p.value) of at least min.probability. marginal.probability is a lower bound on the marginal probability that an alternative hypothesis is true; the default of 0 is conservative.
{
  if(!is.numeric(marginal.probability)) 
  {
  
  
pi1.iterate<- marginal.probability(p.values=p.values, min.probability=0.5, max.iteration=10, tolerance=tolerance)
	 cat('Using marginal.probability estimate of', pi1.iterate, '\n')
	 return(alternative.probable(p.values = p.values, min.probability =min.probability, marginal.probability = pi1.iterate, max.iteration=1 ,tolerance=tolerance, plot.relative.gain = plot.relative.gain))
  }
  if(!is.numeric(p.values)) stop('p.values should be numeric') # vector class not needed
  if(length(min.probability) != 1 || !is.finite(min.probability)) stop('min.probability should be a single number')
  if(length(marginal.probability) != 1 || !is.finite(marginal.probability)) stop('marginal.probability should be a single number')
  is.probability <- function(prob){is.finite(prob) & prob >= 0 & prob <= 1}
 pvals<- na.omit(p.values)   
 if( !all(is.probability(pvals))) stop('Each element of p.values must be between 0 and 1')
  if(!is.probability(min.probability)) stop('min.probability must be between 0 and 1')
  if(min.probability >= 1) stop('min.probability must be less than 1')
  if(!is.probability(marginal.probability)) stop('marginal.probability must be between 0 and 1')
  
  p <- 1 - min.probability
  pi0 <- 1 - marginal.probability
  alphas <- unique(pvals) # significance levels (Type I error rates)
  
  relative.gain <- function(.alpha, .pvals, .pi0, .p)
  {
  if (call.browser)  browser()
  nrej<- sapply(1:length(.alpha), function(j) {sum(.pvals <=.alpha[j])})
    dFDR <- function(.alpha, .pvals, .pi0, .nrej)
	 # estimate of the decisive false discovery rate
    {
      if(length(.alpha) != length(.nrej)) stop('Please report .alpha-.nrej error to www.davidbicke.com')
      dfdr <- ifelse(.nrej > 0, .pi0 * .alpha / (.nrej / length(.pvals)), 0)
      if(any(!is.finite(dfdr) | dfdr < 0)) stop('Please report dFDR error to www.davidbickel.com')
      ifelse(dfdr <= 1, dfdr, 1)
    }

    .dFDR <- dFDR(.pvals = .pvals, .alpha = .alpha, .pi0 = .pi0, .nrej = nrej)
    nrej * (1 - .dFDR / .p) # nrej * (1 - (.pi0 * .alpha * length(.pvals) / nrej) / .p) == nrej - length(.pvals) * .pi0 * .alpha / .p
  }
  gains <- relative.gain(.alpha = alphas, .pvals = pvals, .pi0 = pi0, .p = p)
  if(!all(is.finite(alphas) & is.finite(gains))) {print('Likely error. Type c to continue or Q to quit.'); browser()}
  max.gain <- max(gains)
  if(plot.relative.gain) plot(alphas, gains, xlab = 'test-wise Type I error', ylab = paste('relative gain (max==', max.gain, ')'))
  if(length(alphas) != length(gains)) stop('Please report alphas-gains error to www.davidbickel.com')
  optimal.alphas <- alphas[max.gain == gains & gains > 0]
  optimal.alpha <- if(length(optimal.alphas) > 0) max(optimal.alphas) else -Inf
  if(optimal.alpha >= 0 && optimal.alpha != min(optimal.alphas)) warning('More than one significance level is optimal.')
 k<- rep(logical (length(p.values)))
 k<- sapply(1:length(p.values), function(j){ (p.values[j]<= optimal.alpha)})
 k
 }
 
marginal.probability<- function(p.values, min.probability = 0.5, max.iteration=10, tolerance= get.marginal.probability.tolerance(), verbose=FALSE)
   {
   #(if(!is.numeric(marginal.probability)) 
      circle<- function(pi0,pi1,tolerance)
      {
     	  circle<- abs(pi0-pi1)<tolerance
     	  circle
      }
      pi1<- 0
   	 for (i in 1:max.iteration)
   	 {
   	  pi0=pi1
       if(verbose)
 						        message(" iteration " , i, " of up to " , max.iteration , " started on " ,  date())
   	 ap0 <- alternative.probable(p.values = p.values, min.probability = 0.5,  marginal.probability = pi0, max.iteration=1, tolerance=tolerance, plot.relative.gain = FALSE)
       ap0 <-na.omit(ap0)
       pi1 <- sum(ap0) / length(ap0)
       if(verbose)
 						        message(" iteration " , i, " of up to " , max.iteration , " finished on " ,  date())
        
        if (circle(pi0,pi1,tolerance))  break
            }
            if (!circle(pi0,pi1,tolerance))  warning('Failed to converge, you may need to increase max.iteration.')
            pi1
   
}