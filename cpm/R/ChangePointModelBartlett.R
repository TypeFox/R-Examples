#change point model using the Bartlett statistic to detect shifts in the variance 
#of a (preferably) Gaussian distribution

setClass(
	Class="ChangePointModelBartlett",
		representation=representation(
				useCox='logical'		
		),
		prototype(
		),
		contains='ChangePointModel'
)



makeChangePointModelBartlett <- function(hs=numeric(),startup=20){
	windowStatistic <- list()
	windowStatistic$S <- numeric()
	windowStatistic$W <- numeric()

	return(new(Class='ChangePointModelBartlett',
		windowStatistic=windowStatistic,
		hs=hs,
        changeDetected=FALSE,
        startup=startup
	))
}


setMethod(f='updateWindowStatistic',signature='ChangePointModelBartlett',definition=
function(cpm,x) {
	S <- x; W <- 0; 
	n <- cpm@n

	if(length(cpm@windowStatistic$S) > 0) {
		oldS <- cpm@windowStatistic$S[length(cpm@windowStatistic$S)]
		oldW <- cpm@windowStatistic$W[length(cpm@windowStatistic$W)]
		S <- oldS + x
		W <- oldW + ((n-1)*(x-oldS/(n-1))^2) / n 
	}
	
	cpm@windowStatistic$S <- c(cpm@windowStatistic$S,S)
	cpm@windowStatistic$W <- c(cpm@windowStatistic$W,W)
	return(cpm@windowStatistic)
})



setMethod(f='getTestStatistics',signature='ChangePointModelBartlett',definition=
function(cpm) {
	S <- cpm@windowStatistic$S
	W <- cpm@windowStatistic$W

	if (cpm@n < cpm@startup) {
		results <- list();
		results$val <- -Inf
		results$MLE <- 0
		results$Ds <- numeric()
		return(results)
	}
	
	n <- cpm@n

	Ds <- numeric(length(S)) 		
			
	res <-.C('cpmMLEBartlett',as.double(S),as.double(W),as.integer(n),Ds=as.double(Ds))
	Ds <- res[['Ds']]	
	len <- length(Ds)
	Ds[c(1,len-1,len-2)] <- 0
	
	results <- list()
	results$val <- max(Ds)	 	
	results$MLE <- which.max(Ds)	
	results$Ds <- Ds
	return(results)
})

