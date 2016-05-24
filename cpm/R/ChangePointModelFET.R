#change point model using fishers exact test for detecting
#a change in a binary stream

setClass(
	Class="ChangePointModelFET",
		representation=representation(
				lambda='numeric'
		),
		prototype(
		),
		contains='ChangePointModel'
)



makeChangePointModelFET <- function(hs=numeric(),startup=20,lambda=0.1) {
	windowStatistic <- list()
	windowStatistic$S <- NULL #stores how many 1s have occurred to left of this point (inclusive)
	windowStatistic$N <- NULL #stores how many points have been seen to left of window (inclusive)

	return(new(Class='ChangePointModelFET',
		windowStatistic=windowStatistic,
		lambda = lambda,
		hs=hs,
		startup=startup,
        changeDetected=FALSE
	))
}

last <- function(L) {return(L[length(L)])}
setMethod(f='updateWindowStatistic',signature='ChangePointModelFET',definition=
function(cpm,x) {
	S <- x; 
	N <- 1;

	if(length(cpm@windowStatistic$S) > 0) {
		S <- last(cpm@windowStatistic$S) + x
		N <- last(cpm@windowStatistic$N) + 1
	}

	cpm@windowStatistic$S <- c(cpm@windowStatistic$S,S)
	cpm@windowStatistic$N <- c(cpm@windowStatistic$N,N)
	return(cpm@windowStatistic)
})



setMethod(f='getTestStatistics',signature='ChangePointModelFET',definition=
function(cpm) {
	S <- cpm@windowStatistic$S
	N <- cpm@windowStatistic$N
	n <- cpm@n
	Ds <- rep(1,length(S))		


	if (length(S) < 20) {
		results <- list();
		results$val <- -Inf
		results$MLE <- 0
		results$Ds <- numeric()
		return(results)
	}

	res <-.C('cpmMLEFET',as.double(S),as.integer(length(S)),as.double(N),as.integer(length(N)),as.integer(n),as.double(cpm@lambda),Ds=as.double(Ds))
	Ds <- res[['Ds']]
	len<-length(Ds)
	Ds[c(1,len-1,len)] <- 0

	results <- list()
	results$val <- max(Ds)	 
	results$MLE <- which.max(Ds)
	results$Ds <- Ds
	return(results)
})



