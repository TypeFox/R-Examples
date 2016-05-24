kmweight <-
function(Y, delta)
	{
	srt<-order(Y)
	sy<-as.double(Y[srt])
	sdelta<-as.integer(delta[srt])
	n <- length(sdelta)
	if(n != length(sdelta) || n != length(Y))
	stop("dimensions of Y and delta don't match!")
	kmweights <- numeric(n)
	kmweights[1] <- 1/n
	for(i in 2 : n)
	 {
	 kmweights[i] <- kmweights[i-1] * (n-i+2)/(n-i+1) * (((n-i+1)/(n-i+2))^sdelta[i-1])
	 }
	kmwts<-kmweights*sdelta
	if (sdelta[n]==0)
	kmwts[n]<-1-sum(kmwts)
	return(kmwts)
}
