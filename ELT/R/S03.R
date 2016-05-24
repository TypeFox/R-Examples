

## ------------------------------------------------------------------------ ##
##  Script R for  "Constructing Entity Specific Prospective Mortality Table ##
##                 Adjustment to a reference"                               ##
## ------------------------------------------------------------------------ ##
##  Script        S03.R                                                     ##
## ------------------------------------------------------------------------ ##
##  Description   Execute method 1                                          ##
## ------------------------------------------------------------------------ ##
##  Authors       Tomas Julien, Frederic Planchet and Wassim Youssef        ##
##                julien.tomas@univ-lyon1.fr                                ##
##                frederic.planchet@univ-lyon1.fr                           ##
##                wassim.g.youssef@gmail.com                                ##
## ------------------------------------------------------------------------ ##
##  Version       01 - 2013/11/06                                           ##
## ------------------------------------------------------------------------ ##

## ------------------------------------------------------------------------ ##
##  FctMethod1 function                                                     ##
## ------------------------------------------------------------------------ ##	
## ---------- Supress SMR == 0 and is.na(SMR) == T, take the mean


FctMethod1 = function(d, e, qref, x1, x2, t1, t2){
	SMR <- sum(d[x1 - min(as.numeric(rownames(d))) + 1,]) / sum(e[x1 - min(as.numeric(rownames(e))) + 1,] * qref[x1 - min(x2) + 1, as.character(t1)])
	QxtFitted <- SMR * qref[, as.character(min(t1) : max(t2))]
	colnames(QxtFitted) <- min(t1) : max(t2)
	rownames(QxtFitted) <- x2
	return(list(SMR = SMR, QxtFitted = QxtFitted, NameMethod = "Method1"))
	}

## ------------------------------------------------------------------------ ##
##  Method1 function                                                        ##
## ------------------------------------------------------------------------ ##



Method1 = function(MyData, AgeRange, Plot = F, Color = MyData$Param$Color){
	M1 <- vector("list", length(MyData)-1)
	names(M1) <- names(MyData)[1:(length(MyData)-1)]
	print("Execute method 1 ...")
	print("Compute the SMR ...")
	for (i in 1 : (length(MyData)-1)){
		.WarningInvalidAge(MyData[[i]]$Dxt, MyData[[i]]$Ext, AgeRange, MyData[[i]]$AgeRef, MyData[[i]]$YearCom)
		M1[[i]] <- FctMethod1(MyData[[i]]$Dxt, MyData[[i]]$Ext, MyData[[i]]$QxtRef, AgeRange, MyData[[i]]$AgeRef, MyData[[i]]$YearCom, MyData[[i]]$YearRef)
		print(paste("QxtFitted", names(MyData)[i]," = ",M1[[i]]$SMR," * QxtRef", names(MyData)[i], sep = ""))
		M1[[i]]=c(M1[[i]],list(AgeMethod=AgeRange))
		}
	if(Plot == T){
		.PlotMethod(M1, MyData, AgeRange, Color)
		}
	return(M1)
	}
	
	