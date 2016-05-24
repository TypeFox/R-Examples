shannon <-function(x){
	#Calculate Shannon's Diversity Index given a list of binomial alleles
	t <- table(x)
	if (length(t)<2) return(NaN)
	P=t[1]/sum(t)
	Q=t[2]/sum(t)
	S=-(P*log(P)+Q*log(Q))
	return(S)
}

