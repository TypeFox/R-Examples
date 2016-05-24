setClass(
	Class="ChangePointModelKS",
		representation=representation(
		),
		prototype(
		),
		contains='ChangePointModel'
)

makeChangePointModelKS <- function(hs=numeric(),startup=20) {
	windowStatistic <- list()
	windowStatistic$X <- numeric()
	windowStatistic$N <- numeric()

	return(new(Class='ChangePointModelKS',
        windowStatistic=windowStatistic,
		startup=startup,
        changeDetected=FALSE,
        hs=hs
	))
}

setMethod(f='updateWindowStatistic',signature='ChangePointModelKS',definition=
function(cpm,x) {
    cpm@windowStatistic$X <- c(cpm@windowStatistic$X,x)
	cpm@windowStatistic$N <- c(cpm@windowStatistic$N,cpm@n)
	return(cpm@windowStatistic)
})

setMethod(f='getTestStatistics',signature='ChangePointModelKS',definition=
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
	
	Ds <- rep(Inf,length(X)) 
	len <- length(X)
	
	res<-.C('cpmMLEKS',as.double(X),as.integer(length(X)),as.integer(ord),as.integer(1),as.integer(1),Ds=double(length(X)))
	Ds <- 1-res[['Ds']]
	len <- length(Ds)
	Ds[c(1,len-1,len)] <- 0

	results <- list() 	
	results$val <-  max(Ds) 	 
	results$MLE <- which.max(Ds)
	results$Ds <- Ds
	return(results)
})

