

## ------------------------------------------------------------------------ ##
##  Script R for  "Constructing Entity Specific Prospective Mortality Table ##
##                 Adjustment to a reference"                               ##
## ------------------------------------------------------------------------ ##
##  Script        S10.R                                                     ##
## ------------------------------------------------------------------------ ##
##  Description   Criteria assessing the plausibility and coherence of the  ##
##                mortality trends                                          ##
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
##  Survival Index                                                          ##
## ------------------------------------------------------------------------ ##

.SurvivalFct = function(Qx,l_0,x) {
	Lx <- vector(,length(x))
	Lx[1] <- l_0
	for (i in 1:length(x)) {
		Lx[1+i] <- (1-Qx[i,i])*Lx[i]
	}
	Sx <- Lx / l_0
	return(Sx = Sx[1:length(x)])
}

## ------------------------------------------------------------------------ ##
##  Decomposition                                                           ##
## ------------------------------------------------------------------------ ##

.GetSxx <- function(s, a) {
	x.xx <- (0 : 99) / 100
	nrs <- nrow(s); ncs <- ncol(s)
	a.x <- a[1] : a[2]
	Sxx <- matrix(, 100 * (length(a.x) - 1), nrs)
	a.xx <- vector(, 100 * (length(a.x) - 1))
	for (i in 1 : nrs) {
		l <- 0
		for (j in 1 : (length(a.x) - 1)) {
			for (k in 1 : 100) {
				Sxx[l + k, i] <- (1 - x.xx[k]) * s[i,j] + x.xx[k] * s[i, j + 1] 
				a.xx[l + k] <- a.x[j] + x.xx[k]
			} 
			l <- l + 100
		}
	}
	return(list(Sxx = Sxx, a.xx = a.xx))
}

## ------------------------------------------------------------------------ ##
##  Single indices                                                          ##
## ------------------------------------------------------------------------ ##

.FctSingleIndices = function(q , Age, AgeRef, t1, t2, NameMethod){
	## --- Survival functions
	LSF <- floor(length(min(t1) : max(t2)) / 10) * 10 - 1
	AgeSFMat <- min(AgeRef) : (min(AgeRef) + LSF)
	NAgeSF <- AgeSFMat + 10
	repeat{
		if(max(NAgeSF) > 130){ break }
		AgeSFMat <- rbind(AgeSFMat, NAgeSF)
		NAgeSF <- NAgeSF + 10
		}
	Start  <- AgeSFMat[, 1]
	NberLignes <- nrow(AgeSFMat) - length(Start[Start >= min(Age)])
#	AgeSFMat <- AgeSFMat[(NberLignes + 1) : nrow(AgeSFMat), ]
	Sx <- matrix(, ncol(AgeSFMat), nrow(AgeSFMat))
	for(i in (1+NberLignes) : nrow(AgeSFMat)){
		Sx[,i] <- .SurvivalFct(q[AgeSFMat[i, ] - min(Age) + 1, ], 100, AgeSFMat[i, ])
		}
	## --- Median age at death
	print("Median age at death ...")
	NameMedian <- vector(, nrow(AgeSFMat))
	RangeIndices <- t(apply(AgeSFMat, 1, range))
	for(i in 1 : nrow(AgeSFMat)){
		NameMedian[i] <- paste("Med[", LSF+1, "_T_",RangeIndices[i,1], "]", sep = "")
		}
	Sxx <- vector("list", nrow(AgeSFMat))
	for(i in 1 : nrow(AgeSFMat)){
		Sxx[[i]] <- vector("list",	2)
		Sxx[[i]][[1]] <- Sxx[[i]][[2]] <- vector(, LSF * 100)
		Temp <- .GetSxx(t(Sx[, i]), c(1, LSF))
		Sxx[[i]][[1]] <- Temp$a.xx
		Sxx[[i]][[2]] <- Temp$Sxx
		}
	MedianAge <- matrix(, nrow(AgeSFMat), 1)
	rownames(MedianAge) <- NameMedian
	colnames(MedianAge) <- NameMethod
	for (i in (1+NberLignes) : nrow(AgeSFMat)) {
		MedianAge[i] <- Sxx[[i]][[1]][max(which(Sxx[[i]][[2]] >= 0.5)) + 1]
		}
	print(MedianAge)
	## --- Entropy
	print("Entropy ...")
	NameEntropy <- vector(, nrow(AgeSFMat))
	for(i in 1 : nrow(AgeSFMat)){
		NameEntropy[i] <- paste("H[", LSF + 1, "_T_", RangeIndices[i, 1], "]", sep = "")
		}
	Entropy <- matrix(, nrow(AgeSFMat), 1)
	colnames(Entropy) <- NameMethod
	rownames(Entropy) <- NameEntropy
	for(i in (1+NberLignes) : nrow(AgeSFMat)){
		Entropy[i] <- round(- mean(log(Sx[,i]))/ sum(Sx[,i]), 4)
		}
	print(Entropy)
	## --- Cohort life expectancy for cohort in min(t1)
	print(paste("Cohort life expectancy for cohort in",min(t1),"over",LSF+1,"years ..."))
	NameEspGen <- vector(, nrow(AgeSFMat))
	for(i in 1 : nrow(AgeSFMat)){
		NameEspGen[i] <- paste(LSF+1,"_e_",RangeIndices[i,1],sep="")
		}
	CohortLifeExp <- matrix(, nrow(AgeSFMat), 1)
	for(i in (1+NberLignes) : nrow(AgeSFMat)){
		AgeVec <- AgeSFMat[i, ]
		CohortLifeExp[i] <- sum(cumprod(1 - diag(q[AgeVec + 1 - min(Age), ])))
		}
	colnames(CohortLifeExp) <- NameMethod
	rownames(CohortLifeExp) <- NameEspGen
	print(CohortLifeExp)
	return(list(MedianAge = MedianAge, Entropy = Entropy, CohortLifeExp = CohortLifeExp, NameMethod = NameMethod))
	}

## ------------------------------------------------------------------------ ##
##  Cohort life expectancy // for cohort in min(YearCom) over 5 years       ##
## ------------------------------------------------------------------------ ##

.FctCohortLifeExp5 = function(q, Age, t1, t2, NameMethod){
	AgeExp <- min(Age) : (min(Age)+4)
	AgeMat <- matrix(, length(Age) - length(AgeExp) + 1, length(AgeExp))
	AgeMat[1, ] <- AgeExp
	for (i in 2 : (length(Age) - length(AgeExp) + 1)){ 
		AgeMat[i, ] <- AgeMat[i - 1, ] + 1 
		}
	NameCohortLifeExp5 <- vector(, nrow(AgeMat))
	for(i in 1 : nrow(AgeMat)){
		NameCohortLifeExp5[i] <- paste(5, "_e_", AgeMat[i, 1], sep = "")
		}
	CohortLifeExp5 <- matrix(,length(Age) - length(AgeExp) + 1, length(min(t1) : (max(t2) - 4)))
	colnames(CohortLifeExp5) <- as.character(min(t1) : (max(t2) - 4))
	rownames(CohortLifeExp5) <- NameCohortLifeExp5
	for (i in 1 : (length(min(t1) : (max(t2) - 4)))){
		for (j in 1 : nrow(AgeMat)){
			AgeVec <- AgeMat[j, ] + 1 - min(Age)
			CohortLifeExp5[j, i] <- sum(cumprod(1 - diag(q[AgeVec, i : (i + 4)])))
			}
		}
	return(list(CohortLifeExp5 = CohortLifeExp5, NameMethod = NameMethod))
	}

## ------------------------------------------------------------------------ ##
##  Periodic life expectancies                                              ##
## ------------------------------------------------------------------------ ##

.FctPerLifeExp = function(q, Age, a, t1, t2, NameMethod){
	AgeComp <- min(Age) : pmin(a, max(Age))
	PerLifeExp <- matrix(, length(AgeComp), length(min(t1) : max(t2)))
	colnames(PerLifeExp) <- as.character(min(t1) : max(t2))
	rownames(PerLifeExp) <- AgeComp
	for (i in 1 : (length(min(t1) : max(t2)))){
		for (j in 1 : length(AgeComp)){
			AgeVec <- ((AgeComp[j]) : max(AgeComp)) + 1 - min(Age)
			PerLifeExp[j,i] <- sum(cumprod(1 - q[AgeVec, i]))
			}
		}
	return(list(PerLifeExp = PerLifeExp, NameMethod = NameMethod))
	}

## ------------------------------------------------------------------------ ##
##  ValidationLevel3 function                                               ##
## ------------------------------------------------------------------------ ##



ValidationLevel3 = function(FinalMethod, MyData, Plot = F, Color = MyData$Param$Color, Excel = F){
	AgeFinal <- as.numeric(rownames(FinalMethod[[1]]$QxtFinal))
	print("Validation: Level 3 - Criteria assessing the plausibility and coherence of the mortality trends...")
	SI <- vector("list", length(MyData)-1)
	names(SI) <- names(MyData)[1:(length(MyData)-1)]
	for (i in 1 : (length(MyData)-1)){
		print(paste("Singles indices summarizing the lifetime probability distribution ,",names(MyData)[i]," population ..."))
		SI[[i]] <- .FctSingleIndices(FinalMethod[[i]]$QxtFinal, AgeFinal, MyData[[i]]$AgeRef, MyData[[i]]$YearCom, MyData[[i]]$YearRef, FinalMethod[[i]]$NameMethod)
	}
	if(Excel == T){
		.ExportSingleIndiciesinExcel(SI, MyData)
	}
	print("Computation of the cohorts life expectancies over 5 years ...")
	CohortLifeExp5 <- vector("list", length(MyData)-1)
	names(CohortLifeExp5) <- names(MyData)[1:(length(MyData)-1)]
	for (i in 1 : (length(MyData)-1)){	
		CohortLifeExp5[[i]] <- .FctCohortLifeExp5(FinalMethod[[i]]$QxtFinal, AgeFinal, MyData[[i]]$YearCom, MyData[[i]]$YearRef, FinalMethod[[i]]$NameMethod)
	}
	if(Excel == T){
		.ExportCohortLifeExp5inExcel(CohortLifeExp5)
	}	
	if(Plot == T){
		Path <- "Results/Graphics/Validation"
		.CreateDirectory(paste("/",Path,sep=""))
		print(paste("Create the plot of the cohort life expectancies over 5 years in .../", Path," ...", sep=""))	
		for (i in 1 : (length(MyData)-1)){
			png(filename=paste(Path,"/",FinalMethod[[i]]$NameMethod,"-CohortLifeExp5-",names(MyData)[i],".png", sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
			print(SurfacePlot(as.matrix(CohortLifeExp5[[i]]$CohortLifeExp5),"5_e", paste(FinalMethod[[i]]$NameMethod, "- Cohort Life Exp. over 5 years,",names(MyData)[i],"pop."), c(min(AgeFinal), 130, min(MyData[[i]]$YearCom) ,max(MyData[[i]]$YearRef)), Color))
			dev.off()
		}
		if(length(MyData) == 3){
			print(paste("Create the plot representing the ratio of the cohorts life expectancies over 5 years, ",names(MyData)[1]," / ",names(MyData)[2]," .../", Path," ...", sep=""))
			png(filename=paste(Path,"/",FinalMethod[[1]]$NameMethod,"RatioCohortsLifeExp.png", sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
			print(SurfacePlot(as.matrix(CohortLifeExp5[[1]]$CohortLifeExp5/CohortLifeExp5[[2]]$CohortLifeExp5), paste("Ratio of cohorts life exp.",names(MyData)[1],"/",names(MyData)[2]), paste(FinalMethod[[1]]$NameMethod,"- Ratio of cohorts life exp. over 5 years,",names(MyData)[1],"/",names(MyData)[2]), c(min(AgeFinal), 130, min(MyData[[i]]$YearCom) ,max(MyData[[i]]$YearRef)), Color))
			dev.off()
			print(paste("Create the plot representing the ratio between the fitted probabilities of death ",names(MyData)[1]," / ",names(MyData)[2]," .../", Path," ...", sep=""))
			png(filename=paste(Path,"/",FinalMethod[[1]]$NameMethod,"RatioProbaDeath.png", sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
			print(SurfacePlot(as.matrix(FinalMethod[[1]]$QxtFinal/FinalMethod[[2]]$QxtFinal), paste("Ratio fitted prob. of death",names(MyData)[1],"/",names(MyData)[2]), paste(FinalMethod[[1]]$NameMethod,"Ratio fitted prob. of death",names(MyData)[1],"/",names(MyData)[2]), c(min(AgeFinal), 130, min(MyData[[i]]$YearCom) ,max(MyData[[i]]$YearRef)), Color))
			dev.off()
		}
	}
	print("Computation of the periodic life expectancies ...")
	PerLifeExpFitted <- PerLifeExpObserved <- PerLifeExpComb <- vector("list", length(MyData)-1)
	names(PerLifeExpFitted) <- names(PerLifeExpObserved) <- names(PerLifeExpComb) <- names(MyData)[1:(length(MyData)-1)]
	AgeComp <- min(AgeFinal) : pmin(max(MyData[[1]]$AgeRef), 95)
	YearObs <- range(colnames(MyData[[1]]$Deaths))
	for (i in 1 : (length(MyData)-1)){	
		PerLifeExpFitted[[i]] <- .FctPerLifeExp(FinalMethod[[i]]$QxtFinal, AgeFinal, max(AgeComp), MyData[[i]]$YearCom, MyData[[i]]$YearRef, FinalMethod[[i]]$NameMethod)
		QObserved <- as.matrix((MyData[[i]]$Deaths/MyData[[i]]$Expo)[AgeComp+1, ])	
		PerLifeExpObserved[[i]] <- .FctPerLifeExp(QObserved, AgeFinal, max(AgeComp) , YearObs[1], YearObs[2], FinalMethod[[i]]$NameMethod)$PerLifeExp
		PerLifeExpComb[[i]] <- cbind(PerLifeExpObserved[[i]],PerLifeExpFitted[[i]]$PerLifeExp[,(length(MyData[[i]]$YearCom)+1):length(min(MyData[[i]]$YearCom):max(MyData[[i]]$YearRef))])
		colnames(PerLifeExpComb[[i]]) <- as.character(c(colnames(MyData[[i]]$Deaths),(max(MyData[[i]]$YearCom)+1):max(MyData[[i]]$YearRef)))
		rownames(PerLifeExpComb[[i]]) <- AgeComp		
	}
	if(Excel == T){
		.ExportPeriodicLifeExpinExcel(PerLifeExpFitted, AgeComp)
	}
	if(Plot == T){
		print(paste("Create graphics representing the ratios of the periodic life expectancies (until ", max(AgeComp) ," years) observed / fitted in .../", Path," ...", sep=""))
		for (i in 1 : (length(MyData)-1)){
			png(filename=paste(Path,"/",FinalMethod[[i]]$NameMethod,"RatioPerLifeExp",names(MyData)[i],".png", sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
			if(length(MyData[[i]]$YearCom)>1){
				print(SurfacePlot(as.matrix(PerLifeExpObserved[[i]][,as.character(MyData[[i]]$YearCom)]/(PerLifeExpFitted[[i]]$PerLifeExp[,as.character(MyData[[i]]$YearCom)])),"Ratio of the periodic life expectancies", paste(FinalMethod[[i]]$NameMethod,"Ratio per. life exp. obs./fitted,",names(MyData)[i],"pop."), c(min(AgeComp), max(AgeComp), min(MyData[[i]]$YearCom),max(MyData[[i]]$YearCom)),Color))
			}
			if(length(MyData[[i]]$YearCom)==1){
				print(plot(AgeComp,as.matrix(PerLifeExpObserved[[i]][,as.character(MyData[[i]]$YearCom)]/(PerLifeExpFitted[[i]]$PerLifeExp[,as.character(MyData[[i]]$YearCom)])), ylab="Ratio of the periodic life expectancy", xlab="Age",main=paste(FinalMethod[[i]]$NameMethod,"- Ratio per. life exp.  obs./fitted., Year",MyData[[i]]$YearCom,",",names(MyData)[i],"pop."), type="l",col=Color,lwd=2)) 
			}
			dev.off()
		}
		AgeCoh <- c(55, 70, 85)
		print(paste("Create graphics representing the coherence of the fitted periodic life expectancies (until ", max(AgeComp) ," years) for age ", AgeCoh[AgeCoh <= max(AgeComp) & AgeCoh >= min(AgeComp)]," in .../", Path," ...", sep=""))
		for(k in 1:length(AgeCoh)){
			if(AgeCoh[k] <= max(AgeComp) & AgeCoh[k] >= min(AgeComp)){
				for (i in 1 : (length(MyData)-1)){
					png(filename=paste(Path,"/",FinalMethod[[i]]$NameMethod, "-CoherencePerLifeExp-", AgeCoh[k], "-", names(MyData)[i],".png", sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
					print(.PlotPerExp(PerLifeExpObserved[[i]],PerLifeExpComb[[i]], PerLifeExpFitted[[i]], AgeCoh[k], AgeComp ,c(55,145), Color, min(as.numeric(colnames(MyData[[i]]$Deaths))):max(MyData[[i]]$YearRef), c(min(PerLifeExpComb[[i]][AgeCoh[k]-min(AgeComp)+1,])-.5,max(PerLifeExpComb[[i]][AgeCoh[k]-min(AgeComp)+1,])+.5),names(MyData)[i], FinalMethod[[i]]$NameMethod))
					dev.off()
				}
			}
		}
	}
	return(list(SingleIndices = SI, CohortLifeExp5 = CohortLifeExp5, PerLifeExpFitted = PerLifeExpFitted, PerLifeExpObserved = PerLifeExpObserved, PerLifeExpComb = PerLifeExpComb))
}
