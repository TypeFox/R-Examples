

## ------------------------------------------------------------------------ ##
##  Script R for  "Constructing Entity Specific Prospective Mortality Table ##
##                 Adjustment to a reference"                               ##
## ------------------------------------------------------------------------ ##
##  Script        S09.R                                                     ##
## ------------------------------------------------------------------------ ##
##  Description   Criteria assessing the regularity of the fit              ##
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
##  Computation of the quantities and tests for level 2                     ##
## ------------------------------------------------------------------------ ##

.TestsLevel2 = function(x, y, z, vv, AgeVec, YearVec, NameMethod){
	RUNSTEST <- matrix(, 7, 1)
	SIGNTEST <- matrix(, 6, 1)
colnames(RUNSTEST) <- colnames(SIGNTEST) <- NameMethod
	rownames(RUNSTEST) <- c("Nber of runs", "Signs (-)", "Signs (+)", "Xi (abs)", "Threshold", "Hyp", "p.val")
	rownames(SIGNTEST) <- c("Signs (+)", "Signs (-)", "Xi", "Threshold", "Hyp", "p.val")
	xx <- x[AgeVec - min(as.numeric(rownames(x))) + 1, YearVec]
	yy <-  y[AgeVec - min(as.numeric(rownames(y))) + 1, YearVec]
	zz <- z[AgeVec - min(as.numeric(rownames(z))) + 1, YearVec]
	## --- Runs test		
	val.SIGN <- sign(yy / zz - xx)
	fac.SIGN <- factor(as.matrix(val.SIGN))
	sign.neg <- sum(levels(fac.SIGN)[1] == fac.SIGN)
    sign.pos <- sum(levels(fac.SIGN)[2] == fac.SIGN)
	val.RUNS <- 1 + sum(as.numeric(fac.SIGN[-1] != fac.SIGN[-length(fac.SIGN)]))
	mean.run <- 1 + 2 * sign.neg*sign.pos / (sign.neg+sign.pos)
	var.run <- 2 * sign.neg * sign.pos * (2 * sign.neg * sign.pos - sign.neg - sign.pos) / ((sign.neg + sign.pos)^2 * (sign.neg + sign.pos - 1))
	Threshold.RUNS <- qnorm(1 - vv/2)
	xi.RUNS <- (val.RUNS - mean.run) / sqrt(var.run)
 	p.val.RUNS <- 2*(1-pnorm(abs(xi.RUNS)))
 	if(abs(xi.RUNS) <= Threshold.RUNS){ RUNSTEST[6, 1] <- "H0" }
 	else { RUNSTEST[6,1] <- "H1" }
	RUNSTEST[1, 1] <- val.RUNS
	RUNSTEST[2, 1] <- sign.neg
	RUNSTEST[3, 1] <- sign.pos
	RUNSTEST[4, 1] <- round(abs(xi.RUNS), 4)
	RUNSTEST[5, 1] <- round(Threshold.RUNS, 4);
	RUNSTEST[7, 1] <- round(p.val.RUNS, 4)
	## --- Sign test
	SIGNTEST[1, 1] <- sign.pos
	SIGNTEST[2, 1] <- sign.neg
	xi.SIGN <- (abs(sign.pos - sign.neg) - 1) / sqrt(sign.pos + sign.neg)
	Threshold.SIGN <- qnorm(1 - vv/2)
	p.val.SIGN <- 2 * (1 - pnorm(abs(xi.SIGN)))
	if(xi.SIGN <= Threshold.SIGN){ SIGNTEST[5, 1] <- "H0"}
	else { SIGNTEST[5, 1] <- "H1" }
	SIGNTEST[3, 1] <- round(xi.SIGN, 4)
	SIGNTEST[4, 1] <- round(Threshold.SIGN, 4)
	SIGNTEST[6, 1] <- round(p.val.SIGN, 4)
	RSLTS <- vector("list", 2)
	RSLTS[[1]] <- RUNSTEST; RSLTS[[2]] <- SIGNTEST; RSLTS[[3]] <- NameMethod
	names(RSLTS) <- c("Runs test", "Signs test", "NameMethod")
	return(RSLTS)
	}

## ------------------------------------------------------------------------ ##
##  Get Tests level 2                                                       ##
## ------------------------------------------------------------------------ ##

.GetCritLevel2 = function(OutputMethod, MyData, ValCrit, AgeCrit){
	CritLevel2 <- .TestsLevel2(OutputMethod$QxtFitted, MyData$Dxt, MyData$Ext, ValCrit, AgeCrit, as.character(MyData$YearCom), OutputMethod$NameMethod)
	return(CritLevel2)
	}

## ------------------------------------------------------------------------ ##
##  ValidationLevel2 function                                               ##
## ------------------------------------------------------------------------ ##


ValidationLevel2 = function(OutputMethod, MyData, ValCrit, AgeCrit, Excel = F){
	print("Validation: Level 2 - Criteria assessing the regularity of the fit ...")
	CritLevel2 <- vector("list", length(MyData)-1)
	names(CritLevel2) <- names(MyData)[1:(length(MyData)-1)]
	for (i in 1 : (length(MyData)-1)){
		.WarningInvalidAge(MyData[[i]]$Dxt, MyData[[i]]$Ext, AgeCrit, MyData[[i]]$AgeRef, MyData[[i]]$YearCom)
		print(paste("Tests & quantities for ",names(MyData)[i]," population ..."))
		CritLevel2[[i]] <-.GetCritLevel2(OutputMethod[[i]], MyData[[i]], ValCrit, AgeCrit)
		print(CritLevel2[[i]])
		}
	if(Excel == T){
		.ExportValidationL2inExcel(CritLevel2)
		}
	return(CritLevel2)
	}
