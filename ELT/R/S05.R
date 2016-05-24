

## ------------------------------------------------------------------------ ##
##  Script R for  "Constructing Entity Specific Prospective Mortality Table ##
##                 Adjustment to a reference"                               ##
## ------------------------------------------------------------------------ ##
##  Script        S05.R                                                     ##
## ------------------------------------------------------------------------ ##
##  Description   Execute method 3                                          ##
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
##  FctMethod3 function                                                     ##
## ------------------------------------------------------------------------ ##



FctMethod3 = function(d, e, qref, x1, x2, t1, t2){
	DB <- cbind(expand.grid(x1, t1), c(d[x1 - min(as.numeric(rownames(d))) + 1, ]), c(e[x1 - min(as.numeric(rownames(e))) + 1, ]), c(qref[x1 - min(x2) + 1, as.character(t1)]))
	colnames(DB) <-c("Age", "Year", "D_i", "E_i", "mu_i")
	DimMat <- dim(qref[, as.character(min(t1) : max(t2))])
	if(length(t1)<10){
		PoisMod <- glm(D_i ~ as.numeric(log(mu_i)) + as.numeric(Age), family=poisson, data = data.frame(DB), offset = log(DB[,"E_i"])) 
		QxtFitted <- matrix(exp(as.numeric(coef(PoisMod)[1]) + as.numeric(coef(PoisMod)[2]) * as.numeric(log(qref[, as.character(min(t1) : max(t2))])) + (x2) * as.numeric(coef(PoisMod)[3])), DimMat[1], DimMat[2])
	 }
	 if(length(t1)>=10){
	 	PoisMod <- glm(D_i ~ as.numeric(log(mu_i)) + as.numeric(Age) * as.numeric(Year), family = poisson,  data = data.frame(DB), offset = log(DB[,"E_i"]))
	 	DataGrid <- expand.grid(x2, min(t1) : max(t2)) 
		IntGrid <- matrix(DataGrid[, 1] * DataGrid[, 2], length(x2), length(min(t1) : max(t2)))
	 	QxtFitted <- matrix(exp(as.numeric(coef(PoisMod)[1]) + as.numeric(coef(PoisMod)[2]) * as.numeric(log(qref[, as.character(min(t1) : max(t2))])) + (DataGrid[,1]) * as.numeric(coef(PoisMod)[3]) + (DataGrid[,2]) * as.numeric(coef(PoisMod)[4]) + IntGrid * as.numeric(coef(PoisMod)[5])), DimMat[1], DimMat[2])
	 	}
	colnames(QxtFitted) <- as.character(min(t1) : max(t2))
	rownames(QxtFitted) <- x2
	return(list(PoisMod = PoisMod, QxtFitted = QxtFitted, NameMethod = "Method3"))
	}

## ------------------------------------------------------------------------ ##
##  Method3 function                                                        ##
## ------------------------------------------------------------------------ ##


Method3 = function(MyData, AgeRange, Plot = F, Color = MyData$Param$Color){
	M3 <- vector("list", length(MyData)-1)
	names(M3) <- names(MyData)[1:(length(MyData)-1)]
	print("Execute method 3 ...")
	print("Compute the parameters of the Poisson model ...")
	for (i in 1 : (length(MyData)-1)){
		.WarningInvalidAge(MyData[[i]]$Dxt, MyData[[i]]$Lxt, AgeRange, MyData[[i]]$AgeRef, MyData[[i]]$YearCom)
		M3[[i]] <- FctMethod3(MyData[[i]]$Dxt, MyData[[i]]$Ext, MyData[[i]]$QxtRef, AgeRange, MyData[[i]]$AgeRef, MyData[[i]]$YearCom, MyData[[i]]$YearRef)
		print(paste("Summary of the Poisson model,",names(MyData)[i],"population"))
		print(summary(M3[[i]]$PoisMod))
		M3[[i]]=c(M3[[i]],list(AgeMethod=AgeRange))	
	}
	if(Plot == T){
		.PlotMethod(M3, MyData, AgeRange, Color)
		}
	return(M3)
	}
