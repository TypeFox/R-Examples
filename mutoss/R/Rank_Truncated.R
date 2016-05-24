# 
# 
# Author: FrankKonietschke
###############################################################################


ranktruncated <- function(pValues, K, silent = FALSE){
	
	L <- length(pValues)
	
	if (K > L){
		warn1 <- paste("K must be smaller than L")
		stop(warn1)
	}
	
#-----Compute the test statistic-----#
	
	index <- order(pValues)
	rindex <- order(index)
	spval <- pValues[index]
	
	w <- prod(spval[1:K])
	
	
#---Compute the used pvalues--------#
	
	w.pvalues <- pValues[rindex]
	
	p.used <- data.frame(Position=index[1:K], pValue=spval[1:K])
	
	
	
	
	
#--Compute the function awt as in the paper-----#

awt <- function(w,t,K){
	if (w<=t^K){
		s <- c(0:(K-1))
		num1 <- K*log(t)-log(w)
		num2 <- w * sum(num1^s/factorial(s))
	}
	if (w > t^K) {
		num2<- t^K
	}
	return(num2)
}

# ----- Compute now the exact distribution-----#
	
	fac1 <- choose(L, K+1)*(K+1)
	t <- seq(0.001,0.999,0.001)
	terg <- c()
	
	for (i in 1:length(t)){
		terg[i] <- (1-t[i])^(L-K-1)*awt(w,t[i],K)
}

#--------------Compute the p-Value -------------#
distribution <- fac1*mean(terg)

#-Regarding the error of numerical integration-#
	
p1 <- (distribution>1)
distribution[p1] <- 1

#--------------Prepare the output---------------#

	if (! silent)
	{
		cat("#----Rank Truncated Product of P-Values (Dubridge and Koeleman; 2003)   \n\n")
			
	}
result <- data.frame(Statistic = w, p.Value=distribution)
return(list(Used.pValue=p.used, RTP=result))
}








