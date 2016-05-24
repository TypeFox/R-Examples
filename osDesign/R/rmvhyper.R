rmvhyper <-
function(Mk, m)
{
	##
	K <- length(Mk)
	M <- sum(Mk)
	
	##
	if(m > M)
		mk <- Mk
	else
	{
		mk <- rep(0, K)
		subSample <- table(sample(rep(1:K, Mk), m))
		mk[as.numeric(names(subSample))] <- subSample
	}
	
	##
	return(mk)
}
