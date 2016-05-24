

## ------------------------------------------------------------------------ ##
##  Script R for  "Constructing Entity Specific Prospective Mortality Table ##
##                 Adjustment to a reference"                               ##
## ------------------------------------------------------------------------ ##
##  Script        S02.R                                                     ##
## ------------------------------------------------------------------------ ##
##  Description   Import the reference tables and select the common years   ##
## ------------------------------------------------------------------------ ##
##  Authors       Tomas Julien, Frederic Planchet and Wassim Youssef        ##
##                julien.tomas@univ-lyon1.fr                                ##
##                frederic.planchet@univ-lyon1.fr                           ##
##                wassim.g.youssef@gmail.com                                ##
## ------------------------------------------------------------------------ ##
##  Version       01 - 2013/11/06                                           ##
## ------------------------------------------------------------------------ ##

## ------------------------------------------------------------------------ ##
##  Definition of the functions                                             ##
## ------------------------------------------------------------------------ ##

## ------------------------------------------------------------------------ ##
##  Import the reference tables                                             ##
## ------------------------------------------------------------------------ ##



AddReference = function(History, ReferenceMale = NULL, ReferenceFemale = NULL){
	MyData <- vector("list", length(History))
	names(MyData) <- names(History)
	for (i in 1:(length(MyData)-1)){
		MyData[[i]] <-  vector("list", 10)
		names(MyData[[i]]) <- c("Deaths","Expo","Indi","QxtRef", "AgeRef", "YearRef", "Dxt", "Ext", "Lxt", "YearCom")
		MyData[[i]]$Deaths = as.matrix(History[[i]]$Deaths)
		MyData[[i]]$Expo = as.matrix(History[[i]]$Expo)
		MyData[[i]]$Indi = as.matrix(History[[i]]$Indi)
	}
	if(length(names(History)[names(History) == "Female"]) == 1){
		print("Import the female reference table ...")
		MyData$Female$QxtRef <- as.matrix(ReferenceFemale)
		MyData$Female$AgeRef <- as.numeric(rownames(ReferenceFemale))
		MyData$Female$YearRef <- as.numeric(colnames(ReferenceFemale))
	}
	if(length(names(History)[names(History) == "Male"]) == 1){
		print("Import the male reference table ...")
		MyData$Male$QxtRef <- as.matrix(ReferenceMale)
		MyData$Male$AgeRef <- as.numeric(rownames(ReferenceMale))
		MyData$Male$YearRef <- as.numeric(colnames(ReferenceMale))
	}
	MyData$Param=History$Param
	print("Select the years in common ...")
	YearCom <- pmax(min(MyData[[1]]$YearRef), as.numeric(colnames(History[[1]]$Deaths)[1])) : as.numeric(colnames(History[[1]]$Deaths)[ncol(History[[1]]$Deaths)])		
	print("Obtain the observed statistics for the ages and years in common ...")
	for (i in 1 : (length(MyData)-1)){
		for(j in 1:3){
			MyData[[i]][[6+j]] <- as.matrix(History[[i]][[j]][(min(MyData[[i]]$AgeRef)+1):nrow(History[[i]][[j]]), as.character(YearCom)])
			colnames(MyData[[i]][[6+j]]) <- YearCom
			rownames(MyData[[i]][[6+j]]) <- min(MyData[[i]]$AgeRef):(nrow(History[[i]][[j]])-1)
		}
		MyData[[i]]$YearCom <- YearCom
	}
	return(MyData)	
}
