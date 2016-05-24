

## ------------------------------------------------------------------------ ##
##  Script R for  "Constructing Entity Specific Prospective Mortality Table ##
##                 Adjustment to a reference"                               ##
## ------------------------------------------------------------------------ ##
##  Script        S08.R                                                     ##
## ------------------------------------------------------------------------ ##
##  Description   Criteria assessing the proximity between the observations ##
##                and the model                                             ##
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
##  Computation of the quantities and tests for level 1                     ##
## ------------------------------------------------------------------------ ##

.TestsLevel1 = function(x, y, z, vv, AgeVec, YearVec, NameMethod){
	quantities <- LRTEST <- matrix(, 4, 1)
	SMRTEST <- WMPSRTEST <- matrix(, 5, 1)
	RESID <- matrix(, 2, 1)
	colnames(LRTEST) <- colnames(SMRTEST) <- colnames(RESID) <- colnames(WMPSRTEST) <- colnames(quantities) <- NameMethod
	rownames(LRTEST) <- c("Xi", "Threshold", "Hyp", "p.val")
	rownames(SMRTEST) <- c("SMR", "Xi", "Threshold", "Hyp", "p.val")
	rownames(RESID) <- c("Std. Res.  > 2", "Std. Res.  > 3")	
	rownames(quantities) <- c("Chi2", "R2", "MAPE", "Deviance")
	rownames(WMPSRTEST) <- c("W", "Xi", "Threshold", "Hyp", "p.val")
		xx <- x[AgeVec - min(as.numeric(rownames(x))) + 1, YearVec]
		yy <-  y[AgeVec - min(as.numeric(rownames(y))) + 1, YearVec]
		zz <- z[AgeVec - min(as.numeric(rownames(z))) + 1, YearVec]
	## --- LR test
		val.all <- matrix(, length(AgeVec), length(YearVec))
		colnames(val.all) <- as.character(YearVec)
		val.all[(yy > 0)] <- 2 * yy[(yy > 0)] * log(yy[(yy > 0)] / (xx[(yy > 0)] * zz[(yy > 0)])) - (yy[(yy > 0)] - xx[(yy > 0)] * zz[(yy > 0)] )
		val.all[(yy == 0)] <- 2 * xx[(yy == 0)] * zz[(yy == 0)]
		xi.LR <- sum(val.all)
		Threshold.LR <- qchisq(1 - vv, df = length(xx))
		pval.LR <- 1 - pchisq(xi.LR, df = length(xx))
		LRTEST[1, 1] <- round(xi.LR, 2); LRTEST[2, 1] <- round(Threshold.LR,2) ; 
		if(xi.LR <= Threshold.LR){ LRTEST[3, 1] <- "H0" } else { LRTEST[3, 1] <- "H1" }
		LRTEST[4, 1] <- round(pval.LR, 4) 
	## --- SMR test
		d <- sum(yy); e <- sum(zz * xx)
		val.SMR <- d / e 
 		if(val.SMR >= 1){
 		xi.SMR <- abs(3 * d^(1/2) * (1 - (9 * d)^(- 1) - (d / e)^(- 1 / 3))) }
 		if(val.SMR < 1){
		xi.SMR <- abs(3 * (d + 1)^(1 / 2) * (((d + 1)/e)^(- 1 / 3) + (9 * (d + 1))^(- 1) - 1)) }
		Threshold.SMR <- qnorm(1 - vv)
 		p.val.SMR <- 1 - pnorm(xi.SMR)
 		if(xi.SMR <= Threshold.SMR){ SMRTEST[4, 1] <- "H0" } else { SMRTEST[4,1] <- "H1" }
		SMRTEST[1, 1] <- round(val.SMR, 4); SMRTEST[2, 1] <- round(xi.SMR, 4);
		SMRTEST[3, 1] <- round(Threshold.SMR, 4); SMRTEST[5, 1] <- round(p.val.SMR, 4);
	## --- Wilcoxon Matched-Paris Signed-Ranks test
		tab.temp <- matrix(, length(xx), 3)
		tab.temp[, 1] <- yy / zz - xx
		tab.temp[, 2] <- abs(tab.temp[, 1])
		tab.temp[, 3] <- sign(tab.temp[, 1])
		tab2.temp <- cbind(tab.temp[order(tab.temp[, 2]), ], 1 : length(xx))
		W.pos <- sum(tab2.temp[, 4][tab2.temp[, 3] > 0])
		W.neg <- sum(tab2.temp[, 4][tab2.temp[, 3] < 0])
		WW <- pmax(W.pos, W.neg)
		Xi.WMPSRTEST <- (WW - .5 - length(xx) * (length(xx) + 1) / 4) / sqrt(length(xx) * (length(xx) + 1) * (2 * length(xx) + 1) / 24)
		Threshold.WMPSRTEST <- qnorm(1 - vv / 2)
		p.val.WMPSRTEST <- 2 * (1 - pnorm(abs(Xi.WMPSRTEST)))
		if(Xi.WMPSRTEST <= Threshold.WMPSRTEST){ WMPSRTEST[4, 1] <- "H0" } else { WMPSRTEST[4, 1] <- "H1" }	
		WMPSRTEST[1, 1] <- WW
		WMPSRTEST[2, 1] <- round(Xi.WMPSRTEST, 4)
		WMPSRTEST[3, 1] <- round(Threshold.WMPSRTEST, 4)
		WMPSRTEST[5, 1] <- round(p.val.WMPSRTEST, 4)
	## --- Standardized residuals
		val.RESID <- (xx - yy / zz) / sqrt(xx / zz)
		RESID[1, 1] <- length(val.RESID[(abs(val.RESID)) > 2])
		RESID[2, 1] <- length(val.RESID[(abs(val.RESID)) > 3])
	## --- chi^2, R2, MAPE
		quantities[1, 1] <- round(sum(((yy - zz * xx)^2) / (zz * xx * (1 - xx))), 2)
		quantities[2, 1] <- round(1 - (sum((yy / zz - xx)^2) / sum((yy / zz - mean(yy / zz))^2)), 4)
		quantities[3, 1] <- round(sum(abs((yy / zz - xx)/(yy / zz))[yy > 0])/sum(yy > 0) * 100, 2)
		quantities[4, 1] <- round(2*sum(val.all), 2)
		RSLTS <- vector("list", 5)
		RSLTS[[1]] <- LRTEST; RSLTS[[2]] <- SMRTEST; RSLTS[[3]] <- WMPSRTEST;
		RSLTS[[4]] <- RESID; RSLTS[[5]] <- quantities; RSLTS[[6]] <- NameMethod
		names(RSLTS) <- c("Likelihood ratio test", "SMR test", "Wilcoxon Matched-Pairs Signed-Ranks test", "Standardized residuals", "Quantities", "NameMethod")
		return(RSLTS)
}

## ------------------------------------------------------------------------ ##
##  Get Tests level 1                                                       ##
## ------------------------------------------------------------------------ ##

.GetCritLevel1 = function(OutputMethod, MyData, ValCrit, AgeCrit){
	CritLevel1 <- .TestsLevel1(OutputMethod$QxtFitted, MyData$Dxt, MyData$Ext, ValCrit, AgeCrit, as.character(MyData$YearCom), OutputMethod$NameMethod)
	return(CritLevel1)
	}

## ------------------------------------------------------------------------ ##
##  Compute the deviance                                                    ##
## ------------------------------------------------------------------------ ##

.DevFct = function(x, y, z, AgeVec, YearVec){
		DevMat <- matrix(, length(AgeVec), length(YearVec))
		colnames(DevMat) <- YearVec
		xx <- as.matrix(x)[AgeVec - min(as.numeric(rownames(x))) + 1, YearVec]
		yy <-  as.matrix(y)[AgeVec - min(as.numeric(rownames(y))) + 1, YearVec]
		zz <- as.matrix(z)[AgeVec - min(as.numeric(rownames(z))) + 1, YearVec]
		DevMat[(yy > 0)] <- 2 * (yy[(yy > 0)] * log(yy[(yy > 0)] / (xx[(yy > 0)] * zz[(yy > 0)])) - (yy[(yy > 0)] - xx[(yy > 0)] * zz[(yy > 0)] ))
		DevMat[(yy == 0)] <- 2 * xx[(yy == 0)] * zz[(yy == 0)]
	return(DevMat) }

## ------------------------------------------------------------------------ ##
##  Compute the residuals                                                   ##
## ------------------------------------------------------------------------ ##

.ResFct = function(x, y, z, AgeVec, YearVec, DevMat, NameMethod){
		xx <- as.matrix(x)[AgeVec - min(as.numeric(rownames(x))) + 1, YearVec]
		yy <-  as.matrix(y)[AgeVec - min(as.numeric(rownames(y))) + 1, YearVec]
		zz <- as.matrix(z)[AgeVec - min(as.numeric(rownames(z))) + 1, YearVec]
		RespRes <- as.matrix(yy / zz - xx)
		RespRes[RespRes==-Inf]=NA
		PearRes <- as.matrix((yy - xx * zz) / sqrt(xx * zz))
		PearRes[PearRes==-Inf]=NA
		DevRes <- as.matrix(sign(yy / zz - xx) * sqrt(DevMat))
		DevRes[DevRes==-Inf]=NA
		colnames(RespRes) <- colnames(PearRes) <- colnames(DevRes) <- YearVec
		rownames(RespRes) <- rownames(PearRes) <- rownames(DevRes) <- AgeVec
		Residuals <- list(RespRes, PearRes, DevRes, NameMethod)
		names(Residuals) <- c("Response Residuals", "Pearson Residuals", "Deviance Residuals", "NameMethod")
	return(Residuals) }

## ------------------------------------------------------------------------ ##
##  Nber of deaths and pointwise confidence intervals                       ##
## ------------------------------------------------------------------------ ##
	
.FittedDxtAndConfInt = function(q, e, x1, t1, ValCrit, NameMethod){
	DxtFitted <- as.matrix(q[x1 - min(as.numeric(rownames(q))) + 1, as.character(t1)] * e[x1 - min(as.numeric(rownames(e))) + 1, as.character(t1)])
	DIntUp <- as.matrix(DxtFitted + qnorm(1 - ValCrit / 2) * sqrt(DxtFitted * (1 - q[x1 - min(as.numeric(rownames(q))) + 1, as.character(t1)])))
	DIntLow <- as.matrix(DxtFitted - qnorm(1 - ValCrit / 2) * sqrt(DxtFitted * (1 - q[x1 - min(as.numeric(rownames(q))) + 1, as.character(t1)])))
	DIntLow[DIntLow < 0 ] <- 0
	colnames(DIntUp) <- colnames(DIntLow) <- colnames(DxtFitted) <- as.character(t1)
	rownames(DIntUp) <- rownames(DIntLow) <- rownames(DxtFitted) <- x1
	return(list(DxtFitted = DxtFitted, DIntUp = DIntUp, DIntLow = DIntLow, NameMethod = NameMethod))
	}
	
## ------------------------------------------------------------------------ ##
##  ValidationLevel1 function                                               ##
## ------------------------------------------------------------------------ ##



ValidationLevel1 = function(OutputMethod, MyData, ValCrit, AgeCrit, Plot = F, Color = MyData$Param$Color, Excel = F){
	print("Validation: Level 1 - Criteria assessing the proximity between the observations and the model ...")
	CritLevel1 <- Residuals <- DeathsIntConf <- vector("list", length(MyData)-1)
	names(CritLevel1) <- names(Residuals) <- names(DeathsIntConf) <- names(MyData)[1:(length(MyData)-1)]
	for (i in 1 : (length(MyData)-1)){
		.WarningInvalidAge(MyData[[i]]$Dxt, MyData[[i]]$Ext, AgeCrit, MyData[[i]]$AgeRef, MyData[[i]]$YearCom)
		print(paste("Tests & quantities for ",names(MyData)[i]," population ..."))
		CritLevel1[[i]] <- .GetCritLevel1(OutputMethod[[i]], MyData[[i]], ValCrit, AgeCrit)
		print(CritLevel1[[i]])
		Dev <- .DevFct(OutputMethod[[i]]$QxtFitted, MyData[[i]]$Dxt, MyData[[i]]$Ext, AgeCrit, as.character(MyData[[i]]$YearCom))
		Residuals[[i]] <- .ResFct(OutputMethod[[i]]$QxtFitted, MyData[[i]]$Dxt, MyData[[i]]$Ext, AgeCrit, as.character(MyData[[i]]$YearCom), Dev, OutputMethod[[i]]$NameMethod)
		DeathsIntConf[[i]] <- .FittedDxtAndConfInt(OutputMethod[[i]]$QxtFitted, MyData[[i]]$Ext, AgeCrit, MyData[[i]]$YearCom, ValCrit, OutputMethod[[i]]$NameMethod)
		}
	if(Excel == T){
		.ExportValidationL1inExcel(CritLevel1)
		}
	if(Plot == T){
		Path <- "Results/Graphics/Validation"
		.CreateDirectory(paste("/",Path,sep=""))
		print(paste("Create the plot of the fits for the years in common in .../",Path," ...", sep=""))
		for (i in 1 : length(CritLevel1)){
			for(j in MyData[[i]]$YearCom){
				png(filename=paste(Path,"/",OutputMethod[[i]]$NameMethod,"-Fit-",j,"-",names(MyData)[i],".png",sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
				.PlotFittedYear(OutputMethod[[i]], MyData[[i]], AgeCrit, Color, j, names(MyData)[i])
				dev.off()
				png(filename=paste(Path,"/",OutputMethod[[i]]$NameMethod,"-FitLog-",j,"-",names(MyData)[i],".png",sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
				.PlotFittedYearLog(OutputMethod[[i]], MyData[[i]], AgeCrit, Color, j, names(MyData)[i])
				dev.off()
				}
			}
		print(paste("Create the plot of the residuals for the years in common in .../",Path," ...", sep=""))		
		for (i in 1 : length(CritLevel1)){
			for(j in MyData[[i]]$YearCom){
				png(filename=paste(Path,"/",OutputMethod[[i]]$NameMethod,"-Residuals-",j,"-",names(MyData)[i],".png",sep=""), width  = 3400, height = 1200, res=300, pointsize= 12)
				print(.PlotRes(Residuals[[i]], AgeCrit, Color, j, names(MyData)[i]))
				dev.off()
				}
			}
		print(paste("Create the graphics presenting the pointwise confidence intervals on the  fitted deaths in .../",Path," ...", sep=""))
		for (i in 1 : length(CritLevel1)){
			for(j in MyData[[i]]$YearCom){
				png(filename=paste(Path, "/",OutputMethod[[i]]$NameMethod,"-ConfIntDeaths-",j,"-",names(MyData)[i],".png",sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
				print(.PlotDIntConf(DeathsIntConf[[i]], MyData[[i]], AgeCrit, j , Color, names(MyData)[i]))
				dev.off()			
				}
			}
		}
	return(list(CritLevel1 = CritLevel1, Residuals = Residuals, DeathsIntConf = DeathsIntConf))
	}
