# Name   : GT
# Desc   : Generate a discrete Generation Time distribution given a distribution type
#          and specifications (mean, sd)
# Date   : 2011/11/09
# Author : Boelle, Obadia
###############################################################################


# Function declaration

generation.time <- function#Generation Time distribution
### Create an object of class GT representing a discretized Generation Time distribution.

(type=c("empirical","gamma","weibull","lognormal"),##<< Type of distribution.
val=NULL, ##<< Vector of values used for the empirical distribution, or c(mean, sd) if parametric.
truncate=NULL,##<< Maximum extent of the GT distribution.
step = 1,##<< Time step used in discretization.
first.half=TRUE,##<< First probability computed on half period.
p0=TRUE ##<< Is probability on day 0 0
) 

  
# Code
  
{
##details<< How the GT is discretized may have some impact on the shape of the distribution.
## For example, the distribution may be discretized in intervals of 1 time step starting 
## at time 0, i.e. [0,1), [1,2), and so on. Or it may be discretized as [0,0.5), [0.5, 1.5), ... (the default).
##details<< If the GT is discretized from a given continuous distribution, 
## the expected duration of the Generation Time will be less than the nominal, 
## it will be in better agreement in the second discretization.
##details<< If p0 is TRUE (default) then the generation time distribution 
## is set to 0 for day 0.

# all tests required 
	type=match.arg(type)

#if empirical, check values
	if (type=="empirical") {
    	GT=val
		if (any(GT <0))
			stop("Values in 'val' must be positive")
		if (sum(GT) >1)
			warning("Values will be standardized to sum to 1")
 		if (!is.null(truncate)) {
			if (truncate < length(val)) {
	       warning(paste("Empirical distribution truncated at length ",truncate))
	       GT = GT[1:truncate]
			}
		}
# if parametric
	} else {
		if (length(val)<2 )
			stop("val= c(mean,sd) must be provided for parametric GT")
		mean = val[1]
		sd = val[2]
		if (any(c(mean,sd)<=0))
			stop("'mean' and 'sd' must be positive")
		if (is.null(truncate)) { 
			tmax = ceiling(mean+10*sd); # sufficiently large
		} else {
			tmax=truncate
		}
		if (first.half) {
			t.scale = c(0,0.5+c(0:tmax))
		} else {
			t.scale = c(0:tmax)
		} 
		if (type == "gamma") {
			a = mean*mean/(sd*sd) 
			s = sd*sd / mean						
			GT = diff(pgamma(t.scale,shape=a,scale=s))
		#other cases; weibull, ...		
		} else if (type=="lognormal") {
			meanlog= log(mean^2/sqrt(mean^2+sd^2))
			sdlog=sqrt(2)*sqrt(log(sqrt(mean^2+sd^2)/mean))
			GT = diff(plnorm(t.scale,meanlog=meanlog,sdlog=sdlog))
		} else if (type=="weibull") {
			cv <- sd/(mean)
			if (cv < 1e-06) {
				nu <- cv/(sqrt(trigamma(1)) - cv * digamma(1))
				shape <- 1/nu
				scale <- (mean)/(1 + nu * digamma(1))
			} else {
				aa <- log(cv^2 + 1)
				nu <- 2 * cv/(1 + cv)
				repeat {
					gb <- (lgamma(1 + 2 * nu) - 2 * lgamma(1 + nu) - aa)/(2 * (digamma(1 + 2 * nu) - digamma(1 + nu)))
					nu <- nu - gb
					if (abs(gb) < 1e-12) break
				}
				shape <- 1/nu
				scale <- exp(log(mean) - lgamma(1 + nu))
			}
			GT = diff(pweibull(t.scale,shape=shape,scale=scale))
		}
		if (is.null(truncate)) {
			# truncate when GI distribution >0.9999
			##details<< If no truncation is provided, the distribution will be truncated at 99.99 percent probability.
			GT.cum = cumsum(GT)
      if(length(GT.cum[GT.cum>0.9999])!=0){
  			truncate = (GT.cum > 0.9999)*(1:length(GT.cum))
  			truncate=min(truncate[truncate>0])
  			if (truncate == 0) warning(paste('provide truncate larger than ',mean + 10 * sd))
  			GT = GT[1:truncate]
			}
		}
	}
	if (p0==TRUE) GT[1]=0
	
	time = 0:(length(GT)-1)
	GT = GT/sum(GT)
	mu=sum(GT * time) 
	sigma = sqrt(sum(GT*time^2) - mu^2)
	return(structure(list(GT=GT,time=time,mean=mu,sd=sigma),class="R0.GT"))
  
  ### A list with components:
  ### \item{GT}{The probabilities for each time unit, starting at time 0.}
  ### \item{time}{The time at which probabilities are calculated.}
  ### \item{mean}{The mean of the discretized GT.}
  ### \item{sd}{The standard deviation of the discretized GT.}
}
