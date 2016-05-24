#' Calculate fitted probabilities or quantiles from a (weighted) linear pool 
#' 
#' Calculates a linear pool given a set of elicited judgements in a \code{fit}
#' object. Then calculates required probabilities or quantiles from the pooled
#' cumulative distribution function.
#' 
#' Quantiles are calculate by first calculating the pooled cumulative
#' distribution function at 100 points, and then using linear interpolation to
#' invert the CDF.
#' 
#' @usage plinearpool(fit, x, d = "best", w = 1)
#' qlinearpool(fit, q, d = "best", w = 1)
#' @aliases plinearpool qlinearpool
#' @param fit The output of a \code{fitdist} command.
#' @param x A vector of required cumulative probabilities P(X<=x)
#' @param q A vector of required quantiles
#' @param d The distribution fitted to each expert's probabilities. This must
#' either be the same distribution for each expert, or the best fitting
#' distribution for each expert. Options are \code{"normal"}, \code{"t"},
#' \code{"gamma"}, \code{"lognormal"}, \code{"logt"},\code{"beta"},
#' \code{"best"}.
#' @param w A vector of weights to be used in the weighted linear pool.
#' @return A probability or quantile, calculate from a (weighted) linear pool
#' (arithmetic mean) of the experts' individual fitted probability.
#' @author Jeremy Oakley <j.oakley@@sheffield.ac.uk>
#' @examples
#' \dontrun{
#' # Expert 1 states P(X<30)=0.25, P(X<40)=0.5, P(X<50)=0.75
#' # Expert 2 states P(X<20)=0.25, P(X<25)=0.5, P(X<35)=0.75
#' # Both experts state 0<X<100.
#' 
#' v <- matrix(c(30, 40, 50, 20, 25, 35), 3, 2)
#' p <- c(0.25, 0.5, 0.75)
#' myfit <- fitdist(vals = v, probs = p, lower = 0, upper = 100)
#' 
#' plinearpool(myfit, x=c(20, 50, 80))
#' qlinearpool(myfit, q=c(0.05, 0.5, 0.95))
#' 
#' # give more weight to first expert
#' plinearpool(myfit, x=c(20, 50, 80), w=c(0.7, 0.3)) 
#' 
#' # force the use of gamma distributions for each expert
#' qlinearpool(myfit, q=c(0.05, 0.5, 0.95), d="gamma") 
#' }
#' @export
plinearpool <-
function(fit, x, d = "best", w = 1){
	
  if(min(w)<0 | max(w)<=0){stop("expert weights must be non-negative, and at least one weight must be greater than 0.")}
  
	n.experts <- nrow(fit$vals)
  if(length(w)==1){
    w <- rep(w, n.experts)
  }
  
	px <- matrix(0, length(x), n.experts)
	weight <- matrix(w/sum(w), length(x), n.experts, byrow = T)
	for(i in 1:n.experts){
		px[,i] <- expertprobs(fit, x, d, ex = i)
	}
	
	apply(px * weight, 1, sum)
			
}
