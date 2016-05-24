

## ------------------------------------------------------------------------ ##
##  Script R for  "Constructing Entity Specific Prospective Mortality Table ##
##                 Adjustment to a reference"                               ##
## ------------------------------------------------------------------------ ##
##  Script        S07.R                                                     ##
## ------------------------------------------------------------------------ ##
##  Description   Completion by the Denuit & Goderniaux 2005 method         ##
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
##  Denuit & Goderniaux 2005 method                                         ##
## ------------------------------------------------------------------------ ##

## ---------- Empirical evidence that the annual probability of death follows
## ---------- an exponential quadratic polynomial model at high ages

.CompletionDG2005 = function(q, x, y, RangeStart, RangeCompletion, NameMethod) {

## ---------- For each calendar year
	
	x <- min(x) : pmin(max(x),100)
	RangeStart <- RangeStart - min(x)
	RangeCompletion <- RangeCompletion - min(x)
	CompletionValues <- matrix(, length(y), 3)
	rownames(CompletionValues) <- y
	colnames(CompletionValues) <- c("OptAge", "OptR2", "OptCt")
	CompletionMat <- matrix(, RangeCompletion[2] - RangeCompletion[1] + 1, length(y))
	QxtFinal <- matrix(, RangeCompletion[2] + 1, length(y))
	colnames(QxtFinal) <- y
	rownames(QxtFinal) <- min(x):130
	for(j in 1 : length(y)) {

## ---------- Find the optimal starting age in the age range
## ---------- RangeStart = c(RangeStart[1],RangeStart[2])

		R2Mat <- matrix(,RangeStart[2] - RangeStart[1] + 1, 2)
		colnames(R2Mat) <- c("Age", "R2")
		quadratic.q <- vector("list", RangeStart[2] - RangeStart[1] + 1)
		for (i in 0 : (RangeStart[2] - RangeStart[1])) {
			AgeVec <- ((RangeStart[1] + i) : (max(x) - min(x)))
			R2Mat[i + 1, 1] <- AgeVec[1]
			FitLM <- lm(log(q[AgeVec + 1, j]) ~ AgeVec + I(AgeVec^2))
			R2Mat[i + 1, 2] <- summary(FitLM)$adj.r.squared
			quadratic.q[[i + 1]] <- as.vector(exp(fitted(FitLM)))
			}

## ---------- Extract the R2, OptAge and the corresponding fitted values
		
		if(any(is.na(R2Mat[,2])) == F){
			OptR2 <- max(R2Mat[, 2], na.rm = T)
			OptAge <- as.numeric(R2Mat[which(R2Mat[, 2] == OptR2), 1])
			quadratic.q.opt <- quadratic.q[[which(R2Mat[, 2] == OptR2)]]
		}
		if(any(is.na(R2Mat[,2]))){
			OptR2 <- NA
			OptAge <- AgeVec[1]
			quadratic.q.opt <- quadratic.q[[i + 1]]
		}
		
## ---------- Find OptCt

		OptCt <- coef(lm(log(quadratic.q.opt) ~ I((RangeCompletion[2] - (OptAge : (max(x) - min(x))))^2) - 1))

## ---------- Compute the fitted values CompletionVal for the age range
## ---------- RangeCompletion = c(RangeCompletion[1],RangeCompletion[2])

		CompletionVal <- vector(,RangeCompletion[2] - RangeCompletion[1] + 1)
		for (i in 0 : (RangeCompletion[2] - RangeCompletion[1])) {
			CompletionVal[i + 1] <- exp(OptCt * (RangeCompletion[2] - (RangeCompletion[1] + i))^2)
		}

## ---------- Store the results

		CompletionValues[j, 1] <- OptAge+min(x)
		CompletionValues[j, 2] <- OptR2
		CompletionValues[j, 3] <- OptCt
		CompletionMat[, j] <- CompletionVal
		QxtFinal[, j] <- c(q[1 : RangeCompletion[1], j], CompletionVal)
		
		## ---------- Ajustement autour de l'age optimal choisi par moyenne
		## ---------- geometrique pour k = 5	
	
	k <- 5
	QxtSmooth <- vector(,2*k+1)
	for(i in (RangeCompletion[1]-k):(RangeCompletion[1]+k)){
		QxtSmooth[1+i-(RangeCompletion[1]-k)] <- prod(QxtFinal[(i-k):(i+k),j])^(1/(2*k+1))
		}	
	QxtFinal[(RangeCompletion[1]-k):(RangeCompletion[1]+k),j] <- QxtSmooth
	}
	return(list(CompletionValues = CompletionValues, CompletionMat = CompletionMat, QxtFinal = QxtFinal, NameMethod = NameMethod))
	}

## ------------------------------------------------------------------------ ##
##  CompletionA function                                                    ##
## ------------------------------------------------------------------------ ##




CompletionA = function(OutputMethod, MyData, AgeRangeOptMale, AgeRangeOptFemale, BegAgeCompMale, BegAgeCompFemale, Color = MyData$Param$Color, ShowPlot = T){
	Age <- as.numeric(rownames(OutputMethod[[1]]$QxtFitted))
	ModCompletion <- vector("list", length(MyData)-1)
	names(ModCompletion) <- names(MyData)[1:(length(MyData)-1)]
	for (i in 1 : (length(MyData)-1)){
		print(paste("Completion by the Denuit & Goderniaux 2005 method -",names(MyData)[i],"..."))
		if(names(MyData)[i] == "Female"){
			if(max(Age)<max(AgeRangeOptFemale)){
				stop("The maximum of the age range selected to close the table is greater than the maximum of the age range of the table. Please, select a valid age range.", call. = F)
			}
			if(max(Age)<max(BegAgeCompFemale)){
				stop("The starting age for which the fitted probabilities of the death are replaced by the values obtained from the completion modelmaximum is greater than the maximum of the age range of the table. Please, select a valid starting age.", call. = F)
			}	
			AgeRangeOpt <- AgeRangeOptFemale
			BegAgeComp <- BegAgeCompFemale 
		}
		if(names(MyData)[i] == "Male"){
			if(max(Age)<max(AgeRangeOptMale)){
				stop("The maximum of the age range selected to close the table is greater than the maximum of the age range of the table. Please, select a valid age range.", call. = F)
			}
			if(max(Age)<max(BegAgeCompMale)){
				stop("The starting age for which the fitted probabilities of the death are replaced by the values obtained from the completion modelmaximum is greater than the maximum of the age range of the table. Please, select a valid starting age.", call. = F)
			}	
			AgeRangeOpt <- AgeRangeOptMale
			BegAgeComp <- BegAgeCompMale 
		}
	ModCompletion[[i]] <- c(.CompletionDG2005(OutputMethod[[i]]$QxtFitted, Age, min(MyData[[i]]$YearCom) : max(MyData[[i]]$YearRef), AgeRangeOpt, c(BegAgeComp,130), OutputMethod[[i]]$NameMethod),list(AgeRangeOpt=AgeRangeOpt,BegAgeComp=BegAgeComp))
	if(ShowPlot == T){
		print(paste("Create graphics comparing Before/After the completion for the years in common,",names(MyData)[i] ,"pop. ..."))
			for(j in MyData[[i]]$YearCom){
				dev.new();
				.BeforeAfterCompletion(ModCompletion[[i]], OutputMethod[[i]], MyData[[i]], Age, names(MyData)[i], Color, j)
			}
			print(paste("Create the graphics of the completed surfaces,",names(MyData)[i] ,"pop..."))
				dev.new();
			print(SurfacePlot(as.matrix(ModCompletion[[i]]$QxtFinal),expression(widetilde(q)[xt]), paste("Fitted prob. of death After Completion,", names(MyData)[i], ",pop."), c(min(Age),130,min(MyData[[i]]$YearCom),max(MyData[[i]]$YearRef)),Color))
				dev.new();
			print(SurfacePlot(as.matrix(log(ModCompletion[[i]]$QxtFinal)),expression(paste("log ",widetilde(q)[xt])), paste("Fitted prob. of death After Completion (log),", names(MyData)[i], ",pop."), c(min(Age),130,min(MyData[[i]]$YearCom),max(MyData[[i]]$YearRef)),Color))
		}
	}	
	return(ModCompletion)
}

## ------------------------------------------------------------------------ ##
##  CompletionB function                                                    ##
## ------------------------------------------------------------------------ ##



CompletionB = function(ModCompletion, OutputMethod, MyData, Color = MyData$Param$Color, Plot = F, Excel = F){
	Age <- as.numeric(rownames(OutputMethod[[1]]$QxtFitted))
	if(Plot == T){
		PathA <- "Results/Graphics/Completion"
		.CreateDirectory(paste("/",PathA,sep=""))
		print(paste("Create graphics comparing Before/After the completion for the years in common in .../", PathA," ...", sep=""))
		for (i in 1 : (length(MyData)-1)){
			for(j in MyData[[i]]$YearCom){
				png(filename=paste(PathA,"/",OutputMethod[[i]]$NameMethod,"-BeforeAfterCompletion-",j,"-",names(MyData)[i],".png",sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
				.BeforeAfterCompletion(ModCompletion[[i]], OutputMethod[[i]], MyData[[i]], Age, names(MyData)[i], Color, j)
				dev.off()
			}
		}
		if (length(MyData) == 3){
			for(j in MyData[[i]]$YearCom){
				png(filename=paste(PathA,"/",OutputMethod[[i]]$NameMethod,"-FitAfterCompletion2pops-",j,".png",sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
				.FitPopsAfterCompletionLog(ModCompletion, MyData, Age, Color, j)
				dev.off()
			}
		}
		PathB <- "Results/Graphics/FinalTables"
		.CreateDirectory(paste("/",PathB,sep=""))
		print(paste("Create the graphics of the completed surfaces in .../", PathB," ...", sep=""))
		for (i in 1 : (length(MyData)-1)){
			png(filename=paste(PathB,"/",OutputMethod[[i]]$NameMethod,"-AfterCompletion-ProbaFitted",names(MyData)[i],".png",sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
			print(SurfacePlot(as.matrix(ModCompletion[[i]]$QxtFinal), expression(widetilde(q)[xt]), paste("Fitted prob. of death After Completion,",OutputMethod[[i]]$NameMethod,",",names(MyData)[i], ",pop."), c(min(Age),130,min(MyData[[i]]$YearCom),max(MyData[[i]]$YearRef)),Color))
			dev.off()
			png(filename=paste(PathB,"/",OutputMethod[[i]]$NameMethod,"-AfterCompletion-LogProbaFitted",names(MyData)[i],".png",sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
			print(SurfacePlot(as.matrix(log(ModCompletion[[i]]$QxtFinal)),expression(paste("log ",widetilde(q)[xt])), paste("Fitted prob. of death After Completion (log),",OutputMethod[[i]]$NameMethod,",", names(MyData)[i], ",pop."), c(min(Age),130,min(MyData[[i]]$YearCom),max(MyData[[i]]$YearRef)),Color))
			dev.off()
		}
	}
	if(Excel == T){
		.ExportFinalTablesInExcel(ModCompletion, MyData, Age)
	}
	Final <-  vector("list", length(MyData)-1)
	names(Final) <- names(MyData)[1:(length(MyData)-1)]
	for (i in 1 : (length(MyData)-1)){
		Final[[i]] <- list(QxtFinal = ModCompletion[[i]]$QxtFinal, NameMethod = OutputMethod[[i]]$NameMethod,AgeRange=OutputMethod[[i]]$AgeMethod,AgeRangeOpt=ModCompletion[[i]]$AgeRangeOpt,BegAgeComp=ModCompletion[[i]]$BegAgeComp)
	if (OutputMethod[[1]]$NameMethod=="Method4")
		Final[[i]]=c(Final[[i]],OutputMethod[[1]][c("P.Opt","h.Opt")])
	}
	return(Final)		
}





