

## ------------------------------------------------------------------------ ##
##  Script R for  "Constructing Entity Specific Prospective Mortality Table ##
##				 Adjustment to a reference"							   ##
## ------------------------------------------------------------------------ ##
##  Script		S01.R													 ##
## ------------------------------------------------------------------------ ##
##  Description   Functions to compute the number of deaths, number of	  ##
##				individuals and Exposure								  ##
## ------------------------------------------------------------------------ ##
##  Authors	   Tomas Julien, Frederic Planchet and Wassim Youssef		##
##				julien.tomas@univ-lyon1.fr								##
##				frederic.planchet@univ-lyon1.fr						   ##
##				wassim.g.youssef@gmail.com								##
## ------------------------------------------------------------------------ ##
##  Version	   01 - 2013/11/06										   ##
## ------------------------------------------------------------------------ ##

## ------------------------------------------------------------------------ ##
##  Definition of the functions											 ##
## ------------------------------------------------------------------------ ##

## ------------------------------------------------------------------------ ##
##  Compute the number of deaths over the period of observation			 ##
## ------------------------------------------------------------------------ ##

.ComputeNbDeaths = function(t, DateBegObservation, DateEndObservation){
	DateBeg <- as.Date(DateBegObservation, format = "%Y/%m/%d")
	DateEnd <- as.Date(DateEndObservation, format = "%Y/%m/%d")
	i <- (t$DateOut >= DateBeg) & (t$DateOut <= DateEnd)
	t_obs <- subset(t, i == TRUE)
	VEC <- as.vector(0:120)
	Age <- trunc(t_obs$AgeOut)
	for (x in 0 : 120) {
		VEC[ x + 1 ] <- sum(t_obs$Status[Age == x] == "deceased")
	}
	return(VEC)
}	

.VentileDeaths = function(t, DateBeg, DateEnd){
	YearBeg <- as.integer(format(DateBeg, "%Y"))
	YearEnd <- as.integer(format(DateEnd, "%Y"))
	e <- as.vector(0:120)
	.VentileDeaths <- data.frame(e)
	for (i in YearBeg : YearEnd) {
		if (i == YearBeg){ beg <- DateBeg }
		else { beg <- paste(as.character(i),"/01/01", sep = "") }
		if (i == YearEnd) { end <- DateEnd }
		else { end <- paste(as.character(i), "/12/31",sep = "") }
		e <- .ComputeNbDeaths(t, beg, end)
		.VentileDeaths <- cbind(.VentileDeaths, e)
		names(.VentileDeaths)[i - YearBeg + 2] <- paste(beg, "-", end, sep = "") }
		.VentileDeaths$e=NULL
		return(.VentileDeaths)
	}

## ------------------------------------------------------------------------ ##
##  Compute the exposure over the period of observation					 ##
## ------------------------------------------------------------------------ ##

.ComputeExpo = function(t, DateBegObservation, DateEndObservation, withExpo){
	NbDaysYear <- 365.25
	DateBeg <- as.Date(DateBegObservation, format = "%Y/%m/%d")
	DateEnd <- as.Date(DateEndObservation, format = "%Y/%m/%d")
	observed <- (t$DateIn <= DateEnd) & (t$DateOut >= DateBeg)
	CensuredPeriode <- (t$DateOut >= DateEnd)|((t$DateOut < DateEnd) & (t$Status == "other"))
	t_censured <- subset(t,observed & CensuredPeriode)
	t_out <- subset(t,observed & (!CensuredPeriode))
	Expo <- as.vector(matrix(nrow=121, ncol=1, data=0))
	AgeInExpoC <- pmax(t_censured$AgeIn, difftime(DateBeg, t_censured$DateOfBirth, "", "days") / NbDaysYear)
	AgeOutExpoC <- pmin(t_censured$AgeOut, difftime(DateEnd, t_censured$DateOfBirth, "", "days") / NbDaysYear)
	if (!withExpo){ AgeOutExpoC <- AgeInExpoC + trunc(AgeOutExpoC - AgeInExpoC) + 1 }
	for (x in 0:120){ Expo[x + 1] <- sum(pmax(0, pmin(x + 1, AgeOutExpoC) - pmax(x, AgeInExpoC))) }
	if (nrow(t_out)>0) {
		AgeInExpoS <- pmax(t_out$AgeIn, difftime(DateBeg, t_out$DateOfBirth, "", "days") / NbDaysYear)
		AgeOutExpoS <- pmin(t_out$AgeOut, difftime(DateEnd, t_out$DateOfBirth, "", "days") / NbDaysYear)
                AgeOutExpoS = trunc(AgeOutExpoS )+((AgeOutExpoS -trunc(AgeOutExpoS ))>0)*1
		for (x in 0:120){ Expo[x + 1] <- Expo[x + 1] + sum(pmax(0, pmin(x + 1, AgeOutExpoS) - pmax(x, AgeInExpoS))) } }
		return(Expo)
	}
	
	.VentileExposition=function(t, DateBeg, DateEnd){
		YearBeg <- as.integer(format(DateBeg, "%Y"))
		YearEnd <- as.integer(format(DateEnd, "%Y"))
		e <- as.vector(0 : 120)
		.VentileExposition <- data.frame(e)
		for (i in YearBeg:YearEnd){
			if (i == YearBeg){ beg <- DateBeg }
			else { beg <- paste(as.character(i), "/01/01", sep = "") }
			if (i == YearEnd){ end <- DateEnd }
			else { end <- paste(as.character(i), "/12/31", sep = "") }
			e <- .ComputeExpo(t, beg, end, T)
			.VentileExposition <- cbind(.VentileExposition,e)
			names(.VentileExposition)[i - YearBeg + 2] <- paste(beg, "-", end, sep = "") }
			.VentileExposition$e=NULL
			return(.VentileExposition)
		}
		
## ------------------------------------------------------------------------ ##
##  Compute the number of individuals over the period of observation		##
## ------------------------------------------------------------------------ ##

.VentileIndi = function(t, DateBeg, DateEnd){
	YearBeg <- as.integer(format(DateBeg, "%Y"))
	YearEnd <- as.integer(format(DateEnd, "%Y"))
	e <- as.vector(0 : 120)
	.VentileExposition <- data.frame(e)
	for (i in YearBeg:YearEnd){
		if (i == YearBeg){ beg <- DateBeg }
		else { beg <- paste(as.character(i), "/01/01", sep = "") }
		if (i == YearEnd){ end <- DateEnd }
		else { end <- paste(as.character(i), "/12/31", sep = "") }
		e <- .ComputeExpo(t, beg, end, F)
		.VentileExposition <- cbind(.VentileExposition,e)
		names(.VentileExposition)[i - YearBeg + 2] <- paste(beg, "-", end, sep = "") }
		.VentileExposition$e=NULL
		return(.VentileExposition)
	}

## ------------------------------------------------------------------------ ##
##  .GetPortfolio function												   ##
## ------------------------------------------------------------------------ ##

.GetPortfolio = function(Portfolio, DateBegObs, DateEndObs, DateFormat){
	
## ---------- DateOfBirth, DateIn & DateOut

Portfolio$DateOfBirth  <- as.Date(Portfolio$DateOfBirth, DateFormat)
Portfolio$DateIn  <- as.Date(Portfolio$DateIn, DateFormat)
Portfolio$DateOut  <- as.Date(Portfolio$DateOut, DateFormat)

## ---------- Periode of observation

DateBegObs <- as.Date(DateBegObs, DateFormat)
DateEndObs <- as.Date(DateEndObs, DateFormat)

## ---------- Retreate Status when DateOut > DateEndObs 

Portfolio$Status[Portfolio$DateOut > DateEndObs] <- "other"
Portfolio$DateOut <- pmin(Portfolio$DateOut, DateEndObs, na.rm = T)

## ---------- Compute AgeIn & AgeOut

Portfolio$AgeInDays <- as.numeric(difftime(Portfolio$DateIn, Portfolio$DateOfBirth, "", "days"))
Portfolio$AgeOutDays <- as.numeric(difftime(Portfolio$DateOut, Portfolio$DateOfBirth, "", "days"))
Portfolio$AgeOutDays[is.na(Portfolio$AgeOutDays) == T] <- DateEndObs - Portfolio$DateOfBirth[is.na(Portfolio$AgeOutDays) == T]
Portfolio$AgeIn <- Portfolio$AgeInDays / 365.25
Portfolio$AgeOut <- Portfolio$AgeOutDays / 365.25

## ---------- Supress aberrant ages

## --- AgeIn <= 0

Portfolio <- subset(Portfolio,(Portfolio$AgeIn > 0))

## --- AgeIn >= 100

Portfolio <- subset(Portfolio,(Portfolio$AgeIn < 100))

## --- AgeOut <= 0

Portfolio <- subset(Portfolio,(Portfolio$AgeOut > 0))

## --- AgeIn >= AgeOut

Portfolio <- subset(Portfolio,(Portfolio$AgeOut-Portfolio$AgeIn) > 0)	
## ---------- Segmentation Male / Female

NamesGender <- levels(as.factor(Portfolio$Gender))
MyPortfolio <- vector("list", length(NamesGender))
names(MyPortfolio) <- NamesGender
for(i in 1:length(NamesGender)){
	MyPortfolio[[i]] <- subset(Portfolio, Portfolio$Gender == NamesGender[i])
}

## ---------- Return	

return(list(MyPortfolio = MyPortfolio, DateBegObs = DateBegObs, DateEndObs = DateEndObs))

## ---------- End
}

## ------------------------------------------------------------------------ ##
##  .GetStatistics function												  ##
## ------------------------------------------------------------------------ ##

.GetStatistics = function(MyPortfolio, DateBegObs, DateEndObs){
	Deaths <- .VentileDeaths(MyPortfolio, DateBegObs, DateEndObs)
	Expo <- .VentileExposition(MyPortfolio, DateBegObs, DateEndObs)
	Indi <- .VentileIndi(MyPortfolio, DateBegObs, DateEndObs)
	YearBeg <- as.integer(format(DateBegObs, "%Y"))
	YearEnd <- as.integer(format(DateEndObs, "%Y"))
	colnames(Deaths) <- colnames(Expo) <- colnames(Indi) <- c(YearBeg : YearEnd)
	rownames(Deaths) <- rownames(Expo) <- rownames(Indi) <- 0 : 120
	return(list(Deaths = Deaths, Expo = Expo, Indi = Indi))
}

## ------------------------------------------------------------------------ ##
##  GetObsStatistics function											   ##
## ------------------------------------------------------------------------ ##



.GetHistory = function(MyPortfolio, DateBegObs, DateEndObs, DateFormat){
	History <- .GetPortfolio(MyPortfolio, DateBegObs, DateEndObs, DateFormat)
	MyData <- vector("list", length(History$MyPortfolio)+1)
	names(MyData) <- c(names(History$MyPortfolio),"Param")
	for (i in 1 : length(History$MyPortfolio)){
		print(paste("Compute the number of deaths, number of individuals and exposure for the ", names(History$MyPortfolio)[i]," population ..."))
		MyData[[i]] <- .GetStatistics(History$MyPortfolio[[i]], History$DateBegObs, History$DateEndObs)
		MyData[[i]]$Expo[MyData[[i]]$Expo < MyData[[i]]$Deaths] <- MyData[[i]]$Deaths[MyData[[i]]$Expo < MyData[[i]]$Deaths]
		MyData[[i]]$Indi[MyData[[i]]$Indi < MyData[[i]]$Deaths] <- MyData[[i]]$Deaths[MyData[[i]]$Indi < MyData[[i]]$Deaths]
		MyData[[i]]$Expo[MyData[[i]]$Expo == 0] <- 1 
		MyData[[i]]$Indi[MyData[[i]]$Indi == 0] <- 1 
	}
	return(MyData)
}





 NoCompletion = function(OutputMethod, MyData, Color = MyData$Param$Color, Plot = F, Excel = F){
	Age <- as.numeric(rownames(OutputMethod[[1]]$QxtFitted))
	FinalMethod <- OutputMethod
	for(i in 1:(length(MyData)-1)){
		names(FinalMethod[[i]])[2] <- "QxtFinal"
		names(FinalMethod[[i]])[4] <- "AgeRange"
		FinalMethod[[i]]$QxtFinal[FinalMethod[[i]]$QxtFinal > 1] <- 1
		FinalMethod[[i]]$QxtFinal[is.na(FinalMethod[[i]]$QxtFinal)] <- 1
	}
	if(Plot == T){
		PathB <- "Results/Graphics/FinalTables"
		.CreateDirectory(paste("/",PathB,sep=""))
		print(paste("Create the graphics of the completed surfaces in .../", PathB," ...", sep=""))
		for (i in 1 : (length(MyData)-1)){
			png(filename=paste(PathB,"/", FinalMethod[[i]]$NameMethod,"-WithoutCompletion-ProbaFitted",names(MyData)[i],".png",sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
			print(SurfacePlot(as.matrix(FinalMethod[[i]]$QxtFinal), expression(widetilde(q)[xt]), paste("Fitted prob. of death without completion,", FinalMethod[[i]]$NameMethod,",",names(MyData)[i], ",pop."), c(min(Age),130,min(MyData[[i]]$YearCom),max(MyData[[i]]$YearRef)),Color))
			dev.off()
			png(filename=paste(PathB,"/", FinalMethod[[i]]$NameMethod,"-WithoutCompletion-LogProbaFitted",names(MyData)[i],".png",sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
			print(SurfacePlot(as.matrix(log(FinalMethod[[i]]$QxtFinal)),expression(paste("log ",widetilde(q)[xt])), paste("Fitted prob. of death without completion (log),", FinalMethod[[i]]$NameMethod,",", names(MyData)[i], ",pop."), c(min(Age),130,min(MyData[[i]]$YearCom),max(MyData[[i]]$YearRef)),Color))
			dev.off()
		}
	}
	if(Excel == T){
		.ExportFinalTablesInExcel(FinalMethod, MyData, Age)
	}
	return(FinalMethod)	
}



## ------------------------------------------------------------------------ ##
##  ReadHistory function													##
## ------------------------------------------------------------------------ ##



ReadHistory = function(MyPortfolio, DateBegObs, DateEndObs, DateFormat, Plot = F, Color = "#A4072E", Excel = F){
	MyData <- .GetHistory(MyPortfolio, DateBegObs, DateEndObs, DateFormat)
	if(Plot == T){
		.PlotHistory(MyData, Color)
	}
	if(Excel == T){
		.ExportHistoryInExcel(MyData)
	}
	MyData$Param <- list(Color = Color)
	return(MyData)
}
