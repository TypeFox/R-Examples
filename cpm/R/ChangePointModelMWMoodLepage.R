# Change point model using the MW, Mood, and Lepage statistics
# All 3 CPMs are defined in this file (defining them in separate
# files causes warnings during package checking process due to 
# the order in which R processes files

#change point model using the Mann-Whitney statistic.
setClass(
Class="ChangePointModelMW",
representation=representation(
),
prototype(
),
contains='ChangePointModel'
)

makeChangePointModelMW <- function(hs=numeric(),startup=20) {
	windowStatistic <- list()
	windowStatistic$N <- numeric()
    windowStatistic$X <- numeric()
    
	return(new(Class='ChangePointModelMW',
    windowStatistic=windowStatistic,
    hs=hs,
    startup=startup,
    changeDetected=FALSE
	))
}


setMethod(f='updateWindowStatistic',signature='ChangePointModelMW',definition=
function(cpm,x) {
    cpm@windowStatistic$X <- c(cpm@windowStatistic$X,x)
	cpm@windowStatistic$N <- c(cpm@windowStatistic$N,cpm@n)
	return(cpm@windowStatistic)
})





setMethod(f='getTestStatistics',signature='ChangePointModelMW',definition=
function(cpm) {
	X <- cpm@windowStatistic$X
	N <- cpm@windowStatistic$N
	n <- N[length(N)]
    
	if (length(X) < cpm@startup) {
		results <- list();
		results$val <- -Inf
		results$MLE <- 0
		results$Ds <- numeric()
		return(results)
	}
    
	ranks <- rank(X)
    
	#for each possible split point i, calculate MW statistic
	res <-.C('cpmMLEMW',as.double(X),as.integer(length(X)),as.integer(N),as.integer(length(N)),as.integer(ranks),as.integer(length(ranks)),Ds=double(length(X)))
	Ds <- res[['Ds']]
	len<-length(Ds)	
    Ds[c(1,len-1,len)]<-0
    
    
	results <- list() 	
	results$val <- max(abs(Ds)) 	 	#max value of U over window	
	results$MLE <- which.max(abs(Ds))	#MLE of change point
	results$Ds <- abs(Ds)
	return(results)
})




#Change point model using the Mood statistic for a change in scale parameter.


setClass(
Class="ChangePointModelMood",
representation=representation(
),
prototype(
),
contains='ChangePointModel'
)

makeChangePointModelMood <- function(hs=numeric(),startup=20) {
	windowStatistic <- list()
	windowStatistic$N <- numeric()
    windowStatistic$X <- numeric()
	
    return(new(Class='ChangePointModelMood',
    windowStatistic=windowStatistic,
    startup=20,
    changeDetected=FALSE,
    hs=hs
	))
}

setMethod(f='updateWindowStatistic',signature='ChangePointModelMood',definition=
function(cpm,x) {
    cpm@windowStatistic$X <- c(cpm@windowStatistic$X,x)
	cpm@windowStatistic$N <- c(cpm@windowStatistic$N,cpm@n)
	return(cpm@windowStatistic)
})


setMethod(f='getTestStatistics',signature='ChangePointModelMood',definition=
function(cpm) {
	X <- cpm@windowStatistic$X
	N <- cpm@windowStatistic$N
	n <- N[length(N)]
    
	if (length(X) < cpm@startup) {
		results <- list();
		results$val <- 0
		results$MLE <- 0
		results$Ds <- numeric()
		return(results)
	}
    
	ranks <- rank(X)
    
	res <-.C('cpmMLEMood',as.double(X),as.integer(length(X)),as.integer(N),as.integer(length(N)),as.integer(ranks),as.integer(length(ranks)),Ds=double(length(X)))
	Ds <- res[['Ds']]
	len <- length(Ds)
    Ds[c(1,len-1,len)]<-0
	
	results <- list() 	
	results$val <- max(abs(Ds))			
	results$MLE <- which.max(abs(Ds))	#MLE of change point (two-sided)
	results$Ds <- abs(Ds)
	return(results)
})





setClass(
	Class="ChangePointModelLepage",
		representation=representation(
				cpmMood='ChangePointModelMood',
				cpmMW='ChangePointModelMW'
		),
		prototype(
		),
		contains='ChangePointModel'
)

makeChangePointModelLP <- function(hs=numeric(),startup=20) {
	return(new(Class='ChangePointModelLepage',
        cpmMood=makeChangePointModelMood(numeric(),startup),
		cpmMW=makeChangePointModelMW(numeric(),startup),
		hs=hs,
		startup=startup,
        changeDetected=FALSE
    ))
}


setMethod(f='updateWindowStatistic',signature='ChangePointModelLepage',definition=
function(cpm,x) {
	return(cpm@windowStatistic)
})



cpmLepageProcessObservation <- function(cpm,x) {
	cpm@cpmMood <- processObservation(cpm@cpmMood,x)
	cpm@cpmMW <- processObservation(cpm@cpmMW,x)
	return(cpm)
}


setMethod(f='getTestStatistics',signature='ChangePointModelLepage',definition=
function(cpm) {
	if (cpm@n < cpm@startup) {
			results <- list();
			results$val <- -Inf
			results$MLE <- 0
			results$Ds <- numeric()
			return(results)
		
	}
    
	MWDs <- getTestStatistics(cpm@cpmMW)$Ds
	MoodDs <- getTestStatistics(cpm@cpmMood)$Ds
    
	Ds <- MWDs^2 + MoodDs^2

	results <- list() 	
	results$val <- max(Ds) 	 
	results$MLE <- which.max(Ds)	
	results$Ds <- Ds
	return(results)
})	

cpmResetLepage <- function(cpm) {
    cpm@n <- 0
    cpm@changeDetected <- FALSE
    cpm@cpmMW <- cpmReset(cpm@cpmMW)
    cpm@cpmMood <- cpmReset(cpm@cpmMood)
    return(cpm)
}
	
			