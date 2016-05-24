abf <-
function(theta, V, W, numerator=0, pi1=NA){
	if(length(theta)!=length(V))
		stop("theta and V must have the same length.")
	if(any(V<=0))
		stop("V must be a positive value.")
	if(missing(W)){
		warning("Since W is missing, it is set to (log(2)/qnorm(0.95))^2.")
		W <- (log(2) / qnorm(0.95))^2
	}
	if(!length(W) %in% c(1, length(theta)))
		stop("W must be either of length 1 or have the same length as theta.")
	if(!numerator %in% 0:1)
		stop("numerator must be either 0 or 1.")
	vw <- V + W
	#z2 <- theta * theta / V
	#bf <- sqrt(vw/V) * exp(-z2/2 * W/vw)
	p0 <- dnorm(theta, sd=sqrt(V))
	p1 <- dnorm(theta, sd=sqrt(vw))
	bf <- if(numerator==0) p0/p1 else p1/p0
	#if(numerator==1)
	#	bf <- 1/bf
	if(is.na(pi1))
		return(bf)
	if(!length(pi1) %in% c(1, length(theta)))
		stop("pi1 must either be of length 1 or have the same length as theta.") 
	if(pi1>=1 | pi1<=0)
		stop("pi1 must be between 0 and 1.")
	po <- if(numerator==0) (1-pi1) / pi1 else pi1 / (1-pi1)
	posto <- bf * po
	out <- list(ABF=bf, priorOdds=po, postOdds=posto, BFDP=posto / (1+posto))
	if(numerator==1)
		names(out)[4] <- "PPA"
	out
}

