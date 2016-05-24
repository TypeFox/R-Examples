

## ------------------------------------------------------------------------ ##
##  Script R for  "Constructing Entity Specific Prospective Mortality Table ##
##                 Adjustment to a reference"                               ##
## ------------------------------------------------------------------------ ##
##  Script        S04.R                                                     ##
## ------------------------------------------------------------------------ ##
##  Description   Execute method 2                                          ##
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
##  Logit function                                                          ##
## ------------------------------------------------------------------------ ##

.FctLogit = function(q) { log(q / (1 - q)) }

## ------------------------------------------------------------------------ ##
##  FctMethod2 function                                                     ##
## ------------------------------------------------------------------------ ##



FctMethod2 = function(d, e, qref, x1, x2, t1, t2){
	Qxt <- (d[x1 - min(as.numeric(rownames(d))) + 1, ] / e[x1 - min(as.numeric(rownames(e))) + 1, ])
	LogitQxt <- .FctLogit(Qxt)
	LogitQxt[LogitQxt == -Inf] <- 0
	LogitQxtRef <- .FctLogit(qref[x1 - min(x2) + 1, as.character(t1)])
	LogitQxtRef[LogitQxtRef == -Inf] <- 0
	Distance = function(p){
		LogitQxtFit <- (p[1] + p[2] * LogitQxtRef)
		QxtFit <- exp(LogitQxtFit) / (1 + exp(LogitQxtFit))
		sum(abs(e[x1 - min(as.numeric(rownames(e))) + 1, ] * (Qxt - QxtFit)))
		}
	ModPar <- constrOptim(c(0, 1), Distance, ui = c(0, 1), ci = 0, control = list(maxit = 10^3), method = "Nelder-Mead")$par
	LogitQxtRef <-  .FctLogit(qref[, as.character(min(t1) : max(t2))])
	QxtFitted <- as.matrix(exp(ModPar[1] + ModPar[2] * LogitQxtRef) / (1 + exp(ModPar[1] + ModPar[2] * LogitQxtRef)))
	colnames(QxtFitted) <- as.character(min(t1) : max(t2))
	rownames(QxtFitted) <- x2
	return(list(ModPar = ModPar, QxtFitted = QxtFitted, NameMethod = "Method2"))
}

## ------------------------------------------------------------------------ ##
##  Method2 function                                                        ##
## ------------------------------------------------------------------------ ##


Method2 = function(MyData, AgeRange, Plot = F, Color = MyData$Param$Color){
	
	M2 <- vector("list", length(MyData)-1)
	names(M2) <- names(MyData)[1:(length(MyData)-1)]
	print("Execute method 2 ...")
	print("Compute the parameters of the semi-parametric relational model ...")
	for (i in 1 : (length(MyData)-1)){
		.WarningInvalidAge(MyData[[i]]$Dxt, MyData[[i]]$Ext, AgeRange, MyData[[i]]$AgeRef, MyData[[i]]$YearCom)
		M2[[i]] <- FctMethod2(MyData[[i]]$Dxt, MyData[[i]]$Ext, MyData[[i]]$QxtRef, AgeRange, MyData[[i]]$AgeRef, MyData[[i]]$YearCom, MyData[[i]]$YearRef)
		print(paste("logit (QxtFitted", names(MyData)[i],") = ",M2[[i]]$ModPar[1]," + ",M2[[i]]$ModPar[2]," * logit (QxtRef", names(MyData)[i],")", sep = ""))
		M2[[i]]=c(M2[[i]],list(AgeMethod=AgeRange))	
	}
	if(Plot == T){
		.PlotMethod(M2, MyData, AgeRange, Color)
	}
	return(M2)
}
