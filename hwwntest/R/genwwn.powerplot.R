genwwn.powerplot <-
function(N=c(32,64,128,256,512,1024), ar=NULL, ma=NULL,
plot.it=TRUE, sigsq=1, alpha=0.05, away.from="standard", filter.number=10,
family="DaubExPhase", verbose=FALSE, ylim=c(0,1)){

lN <- length(N)

if (lN < 2)
	stop("Vector N has to contain at least two integers.")

answer <- rep(0, lN)

#
# Compute theoretical power approximation for this model for
# the different sample sizes
#

for(i in 1:lN)	{

	answer[i] <- genwwn.thpower(N = N[i], ar = ar, ma = ma,
		sigsq = sigsq, alpha = alpha,
		away.from = away.from, filter.number = filter.number, 
	        family = family, verbose = verbose)$th.power
	}

alist <- list(N=N, power=answer, ar=ar, ma=ma, sigsq=sigsq,
		alpha=alpha, away.from=away.from, filter.number=filter.number,
		family=family)
if (plot.it==TRUE)	{
	plot(N, answer, type="l", xlab="Sample Size",
		ylab="Approximate Power of the genwwn test", ylim=ylim,
		lty=2)
	points(N, answer, pch="x")
	}


return(alist)


}
