#change point model for joint gaussian mean/variance monitoring

setClass(
	Class="ChangePointModelJoint",
		representation=representation(	
		),
		prototype(
		),
		contains='ChangePointModelStudent'
)


makeChangePointModelJoint <- function(hs=numeric(), startup=20) {
	windowStatistic <- list()
	windowStatistic$S <- numeric()
	windowStatistic$W <- numeric()

	return(new(Class='ChangePointModelJoint',
		windowStatistic=windowStatistic,
		hs=hs,
        startup=startup,
        changeDetected=FALSE
	))
}



setMethod(f='getTestStatistics',signature='ChangePointModelJoint',definition=
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
		
	res <- .C('cpmMLEJoint',as.double(S),as.integer(length(S)),as.double(W),as.integer(length(W)),as.integer(n),as.integer(0),Ds=as.double(Ds))
	Ds <- abs(res[['Ds']])
	len<-length(Ds)
	
    Ds[c(1,len-1,len)] <- 0	
	results <- list()   
	results$val <- max(Ds)          #max value of T
	results$MLE <- which.max(Ds)	#MLE of change point
	results$Ds <- abs(Ds)
	return(results)
})

