kmweight.corr <-
function(Y, delta)
	{
	srt<-order(Y)
	sy<-as.double(Y[srt])
	sdelta<-as.integer(delta[srt])
	n <- length(sdelta)
	if(n != length(sdelta) || n != length(Y))
	stop("dimensions of Y and delta don't match!")
	biasf <- numeric(n-2)
      biasf[1]<-((n-2)/(n-1))^sdelta[1]
      for(i in 2 : (n-2))
      {
      biasf[i]<- biasf[i-1]*(((n-i-1)/(n-i))^sdelta[i])
      }
      kmweights <- numeric(n)
	 kmweights[1] <- 1/n
	 for(i in 2 : n)
	 {
	 kmweights[i] <- kmweights[i-1] * (n-i+2)/(n-i+1) * (((n-i+1)/(n-i+2))^sdelta[i-1])
	 }
	 km.n<-kmweights[n]
       kmweights[n] <-km.n+(1-(1/n))*biasf[n-2]
	 kmwts<-kmweights*sdelta
	 return(kmwts)
	}
