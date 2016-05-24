#change point model using the student T statistic to detect shifts in the mean 
#of a (preferably) Gaussian distribution

setClass(
	Class="ChangePointModelStudent",
		representation=representation(	
		),
		prototype(
		),
		contains='ChangePointModel'
)



makeChangePointModelStudent <- function(hs=numeric(), startup=20) {
	windowStatistic <- list()
	windowStatistic$S <- numeric()
	windowStatistic$W <- numeric()

	return(new(Class='ChangePointModelStudent',
		windowStatistic=windowStatistic,
		hs=hs,
        startup=startup,
        changeDetected=FALSE
	))
}


setMethod(f='updateWindowStatistic',signature='ChangePointModelStudent',definition=
function(cpm,x) {
	S <- x; W <- 0; 
	n <- cpm@n

	if(length(cpm@windowStatistic$S) > 0) {
		oldS <- cpm@windowStatistic$S[length(cpm@windowStatistic$S)]
		oldW <- cpm@windowStatistic$W[length(cpm@windowStatistic$W)]
		S <- oldS + x
		W <- oldW + (( (n-1)*x - oldS)^2) /((n*(n-1)))
	}
	
	cpm@windowStatistic$S <- c(cpm@windowStatistic$S,S)
	cpm@windowStatistic$W <- c(cpm@windowStatistic$W,W)
	return(cpm@windowStatistic)
})




setMethod(f='getTestStatistics',signature='ChangePointModelStudent',definition=
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
		
	res <- .C('cpmMLEStudent',as.double(S),as.integer(length(S)),as.double(W),as.integer(length(W)),as.integer(n),as.integer(0),Ds=as.double(Ds))
	Ds <- abs(res[['Ds']])
	len<-length(Ds)
	
    Ds[c(1,len-1,len)] <- 0	
	results <- list()   
	results$val <- max(Ds)          #max value of T
	results$MLE <- which.max(Ds)	#MLE of change point
	results$Ds <- abs(Ds)
	return(results)
})

