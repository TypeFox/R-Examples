fitEcm <-
function(data, max.pop){
	
	data <- classToInd(data)
	
	data <- data[which(data[,2] != 0),]
	
	counts <- data[,2]
	
	cap.ind <- sum(counts)/length(counts)
	sample.size <- sum(counts)
	sampled.ind <- length(counts)
	model <- "Equal.capture"
	
	
## likelihood function for ECM model
	
	likEcm <- function(counts, n){
		
		t.ind <- length(counts)
		
		U <- sum( sapply (1:n, function(u) log(u)))
		
		V <- sum( sapply (1:(n - t.ind), function(v) log(v)))
		
		Y <- (log((1/n)) * (sum(counts)))
		
		if (n == t.ind){
			
			return(U + Y)
			
		} else {
			
			return(U - V + Y)
		}
		
	}

## function to get likelihood constant
	
	getLikConEcm <- function(counts, sample.size){
		c <- - sum(sapply(c(1:length(counts)), function(x) log(x)))

		return(c)
		
	}
	
## if data is all singletons, return ml for n as max.pop
		
	if (all(unique(data[,2]) == 1)){
		warning("Sorry but your data is not informative")
		ml.pop.size <- max.pop
		likelihood <- likEcm(counts, n=max.pop) + getLikConEcm(counts, sample.size)
		return(list(model=model, likelihood=likelihood, ml.pop.size=ml.pop.size, cap.ind=cap.ind, sampled.ind=sampled.ind, sample.size=sample.size, max.pop=max.pop))

	}
	
	vv <- c(length(counts):max.pop)
	
	lik <- sapply(vv, function(x) likEcm(counts, x))
				  
	likelihood <- max(lik) + getLikConEcm(counts, sample.size)
				  
	ml.pop.size	<- vv[which(lik == max(lik))]		  
			
	return(list(model=model, likelihood=likelihood, ml.pop.size=ml.pop.size, cap.ind=cap.ind, sampled.ind=sampled.ind, sample.size=sample.size, max.pop=max.pop))	

}
