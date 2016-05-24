###########################################
# converts probabilities to logits
probs2logits <- function(probs){
	include_zero <- TRUE
	#*** initial values
	NP <- length(probs)
	q1 <- stats::qlogis(probs)
	x0 <- q1-q1[1]
	x0 <- x0[-1]
	#*** define optimization function
	fr <- function(x) {  
		NX <- length(x)  
		NP <- NX+1
		y <- rep(0,NP)      
		y[seq(2,NX+1)] <- x  
		p1 <- exp(y) / sum( exp(y) )
		p1 <- sum( (probs - p1)^2 )
		return(p1)
	}
	#*** optimization
	a1 <- stats::optim( x0 , fr)
	y0 <- a1$par
	if (include_zero){
			y0 <- c(0,y0)
					}
	return(y0)
			}
##################################################			
logits2probs <- function(y){
	x <- exp(y)
	x <- x / sum(x)
	return(x)
		}
##################################################
		