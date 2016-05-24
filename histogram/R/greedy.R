histgreedy<-function(BL,NL,n,binmax,verbose=verbose){
  nB<-length(BL)
	if (nB==2) return(BL)

  #compute possible increses of likelihood by splitting bins

	inclikelihood <- function(i,j) {
		if ((i+1)==j) return(NULL)
		else if (BL[i]==BL[j]) return(rep(-Inf,j-i-1))
		else {		
			k <- (i+1):(j-1)
			old <- (NL[j]-NL[i])*log((NL[j]-NL[i])/(BL[j]-BL[i])/n) 
			new<-rep(0,length(k))
			indv<-(NL[j]-NL[k]>0)
			new[indv]<-((NL[j]-NL[k])*log((NL[j]-NL[k])/(BL[j]-BL[k])/n))[indv]
			indv<-(NL[k]-NL[i]>0)
			new[indv]<-new[indv]+((NL[k]-NL[i])*log((NL[k]-NL[i])/(BL[k]-BL[i])/n))[indv]-old 
			new
			}
		}		
   	breaks <- c(1,nB)
  # compute increment of likelihood obtained by adding one bin
	increment <- c(-Inf,inclikelihood(breaks[1],breaks[2]),-Inf)

	while ((max(increment)>0)&&(length(breaks)<(binmax+1))) {
		maxi <- which.max(increment)
		here <- sum(breaks<maxi)
		i <- breaks[here];
		j <- breaks[here+1];
		debut<-increment[1:i]; debut[i]=-Inf
		fin<-increment[j:nB]; fin[1]=-Inf
		gche <- inclikelihood(i,maxi)
		drte <- inclikelihood(maxi,j)
		increment <- c(debut,gche,-Inf,drte,fin)
		breaks <- c(breaks[1:here],maxi,breaks[(here+1):length(breaks)])
		}
		
		breaks <- BL[breaks]
		if (verbose) message(paste('- Pre-selected finest partition with',length(breaks)-1,'bins.'))
		return(breaks)
}



