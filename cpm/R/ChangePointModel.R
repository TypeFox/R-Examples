#S4 class for change point model object
setClass(
	Class="ChangePointModel",
		representation=representation(
            windowStatistic='list',
			hs='numeric',
			n='numeric',
			startup='numeric',
            changeDetected='logical'
		),
		prototype(
			n=0,
			startup=20,
			changeDetected=FALSE,
			hs=numeric()
		)
)

#returns D_{k,t} statistics
getTestStatistics <- function(cpm) {
	standardGeneric("getTestStatistics")	
}

#updates the CPM given an observation
updateWindowStatistic <- function(cpm,x) {
	standardGeneric("updateWindowStatistic")	
}


setGeneric(name='getTestStatistics',def=getTestStatistics)
setGeneric(name='updateWindowStatistic',def=updateWindowStatistic)