
#change point model using the Cramer-Von-Mises statistic. 

setClass(
	Class="ChangePointModelCVM",
		representation=representation(
		),
		prototype(
		),
		contains='ChangePointModel'
)

makeChangePointModelCVM <- function(hs=numeric(),startup=20) {
	windowStatistic <- list()
    windowStatistic$X <- numeric()
	windowStatistic$N <- numeric()

	return(new(Class='ChangePointModelCVM',
		windowStatistic=windowStatistic,
		hs=hs,
		changeDetected=FALSE,
        startup=startup
	))
}


setMethod(f='updateWindowStatistic',signature='ChangePointModelCVM',definition=
function(cpm,x) {
    cpm@windowStatistic$X <- c(cpm@windowStatistic$X,x)
	cpm@windowStatistic$N <- c(cpm@windowStatistic$N,cpm@n)
	return(cpm@windowStatistic)
})


setMethod(f='getTestStatistics',signature='ChangePointModelCVM',definition=
function(cpm) {
	X <- cpm@windowStatistic$X

	if (length(X) < cpm@startup) {
		results <- list();
		results$val <- -Inf
		results$MLE <- 0
		results$Ds <- numeric()
		return(results)
	}
	ord <- order(X)
    
	len <- length(X)
    res<-.C('cpmMLECVM',as.double(X),as.integer(length(X)),as.integer(ord),Ds=double(length(X)))
	Ds <- res[['Ds']]
    
    Ds[c(1,len-1,len)]<-0

	

	results <- list() 	
	results$val <- max(Ds) 	 	
	results$MLE <- which.max(Ds)
	results$Ds <- Ds
	return(results)
})
