

## ------------------------------------------------------------------------ ##
##  Script R for  "Constructing Entity Specific Prospective Mortality Table ##
##                 Adjustment to a reference"                               ##
## ------------------------------------------------------------------------ ##
##  Script        S06.R                                                     ##
## ------------------------------------------------------------------------ ##
##  Description   Execute method 4                                          ##
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
##  Function for the choice of the optimal parameters                       ##
## ------------------------------------------------------------------------ ##

.GetCrit = function(z, x, hh, pp, k, fam, lk, ww) { 
	AICMAT <- V2MAT <- SpanMat <- matrix(,length(hh),length(pp))
	colnames(AICMAT) <- colnames(V2MAT) <- colnames(SpanMat) <- c("Degree 0", "Degree 1", "Degree 2", "Degree 3")
	rownames(AICMAT) <- rownames(V2MAT) <- rownames(SpanMat) <- c(1:20)
	for ( h in hh ) {
		for( p in (pp+1) ) { 
			Span <- (2*h+1)/length(z)
			SpanMat[h, p] <- Span
			AICMAT[h, p] <- aic(c(z) ~ lp(c(x),deg = p-1, nn = Span,  scale=1), ev=dat(), family = fam, link = lk, kern = k, weights = ww)[4]
			V2MAT[h, p] <-  locfit(c(z) ~ lp(c(x),deg = p-1, nn = Span,  scale=1), ev=dat(), family = fam, link = lk, kern = k, weights = ww)$dp["df2"] } } 
	return(list(AICMAT = AICMAT, V2MAT = V2MAT, SpanMat = SpanMat )) }

## ------------------------------------------------------------------------ ##
##  FctMethod4_1stPart function                                             ##
## ------------------------------------------------------------------------ ##



FctMethod4_1stPart = function(d, e, qref, x1, x2, t1){
	Dx <- apply(as.matrix(d[x1 - min(as.numeric(rownames(d)))+1, ]), 1, sum)
	DxRef <- apply(as.matrix(as.matrix(e)[x1 - min(as.numeric(rownames(e))) + 1, ] * qref[x1 - min(x2) + 1, as.character(t1)]), 1, sum)
	ModCrit <- .GetCrit(Dx, x1, 1 : 20, 0 : 3, c("epan"), "poisson", "log", c(DxRef))
	return(ModCrit)
}

## ------------------------------------------------------------------------ ##
##  Method4A function                                                       ##
## ------------------------------------------------------------------------ ##



Method4A = function(MyData, AgeRange, AgeCrit, ShowPlot = F){
	if(min(AgeRange)>min(AgeCrit)){
		stop("The age range selected to fit Method4 should at least cover the age range used for the validation. Please, select a valid age range.", call. = F)
	}
	if(max(AgeRange)<max(AgeCrit)){
		stop("The age range selected to fit Method4 should at least cover the age range used for the validation. Please, select a valid age range.", call. = F)
	}
	M4A <- vector("list", length(MyData)-1)
	names(M4A) <- names(MyData)[1:(length(MyData)-1)]
	print("Execute method 4 First Part ...")
	print("Choice of the optimal smoothing parameters ...")
	for (i in 1 : (length(MyData)-1)){
		.WarningInvalidAge(MyData[[i]]$Dxt, MyData[[i]]$Lxt, AgeRange, MyData[[i]]$AgeRef, MyData[[i]]$YearCom)
		M4A[[i]] <- FctMethod4_1stPart(MyData[[i]]$Dxt, MyData[[i]]$Ext, MyData[[i]]$QxtRef, AgeRange, MyData[[i]]$AgeRef, MyData[[i]]$YearCom)
		if (ShowPlot == T){
			dev.new();print(.PlotCrit(M4A[[i]]$AICMAT, M4A[[i]]$V2MAT, c(min(M4A[[i]]$AICMAT[3:20,]),max(M4A[[i]]$AICMAT)), c(min(M4A[[i]]$V2MAT[3:20,]),max(M4A[[i]]$V2MAT))))	
			}
		M4A[[i]]=c(M4A[[i]],list(AgeMethod=AgeRange,AgeCrit=AgeCrit))
	}
	return(M4A)
}

## ------------------------------------------------------------------------ ##
##  FctMethod4_2ndPart function                                             ##
## ------------------------------------------------------------------------ ##



FctMethod4_2ndPart = function(d, e, qref, x1, x2, t1, t2, P.Opt, h.Opt){
	Dx <- apply(as.matrix(d[x1 - min(as.numeric(rownames(d)))+1, ]), 1, sum)
	DxRef <- apply(as.matrix(as.matrix(e[x1 - min(as.numeric(rownames(e))) + 1, ]) * qref[x1 - min(x2) + 1, as.character(t1)]), 1, sum)
	QxtFitted <- predict(locfit(Dx ~ lp(x1, deg = P.Opt, nn = (h.Opt * 2 + 1) / length(x1), scale=1), ev = dat(), family = "poisson", link = "log", kern = c("epan"), weights = c(DxRef))) * qref[x1 - min(x2) + 1, as.character(min(t1):max(t2))]
	colnames(QxtFitted) <- as.character(min(t1):max(t2))
	rownames(QxtFitted) <- x1
	return(list(QxtFitted = QxtFitted, NameMethod = "Method4"))
}

## ------------------------------------------------------------------------ ##
##  Method4B function                                                       ##
## ------------------------------------------------------------------------ ##



Method4B = function(PartOne, MyData, OptMale, OptFemale, Plot = F, ShowPlot = F, Color = MyData$Param$Color){
	AgeRange = PartOne[[1]]$AgeMethod
	if(min(AgeRange)>min(PartOne[[1]]$AgeCrit) | max(AgeRange)<max(PartOne[[1]]$AgeCrit)){
		stop("The age range selected to fit Method4 should at least cover the age range used for the validation. Please, select a valid age range.", call. = F)
	}
	M4B <- vector("list", length(MyData)-1)
	names(M4B) <- names(MyData)[1:(length(MyData)-1)]
	for (i in 1 : (length(MyData)-1)){
		print(paste("Execute method 4 Second Part -",names(MyData)[i],"..."))
		if(names(MyData)[i] == "Female"){
			if(OptFemale[1]>3 | OptFemale[1]<0){
				stop("The optimal polynomial degree should be selected between degrees 0 to 3. Please, select a valid optimal polynomial degree.", call. = F)
			}
			if(OptFemale[2]>20 | OptFemale[2]<1){
				stop("The optimal window width should be selected between 1 to 20, i.e. Lambda = value * 2 + 1 observations. Only 3 < Lambda < 41 is allowed. Please, select a window width.", call. = F)
			}
			P.Opt <- OptFemale[1]
			h.Opt <- OptFemale[2]
		}
		if(names(MyData)[i] == "Male"){
			if(OptMale[1]>3 | OptMale[1]<0){
				stop("The optimal polynomial degree should be selected between degrees 0 to 3. Please, select a valid optimal polynomial degree.", call. = F)
			}
			if(OptMale[2]>20 | OptMale[2]<1){
				stop("The optimal window width should be selected between 1 to 20, i.e. Lambda = value * 2 + 1 observations. Only 3 < Lambda < 41 is allowed. Please, select a window width.", call. = F)
			}
			P.Opt <- OptMale[1]
			h.Opt <- OptMale[2] 
		}		
		print("Non-parametric smoothing of the periodic table ...")
		if(ShowPlot == T){		
			dev.new(); .PlotCritChoice(PartOne[[i]]$AICMAT, PartOne[[i]]$V2MAT, c(min(PartOne[[i]]$AICMAT[3:20,]),max(PartOne[[i]]$AICMAT)), c(min(PartOne[[i]]$V2MAT[3:20,]),max(PartOne[[i]]$V2MAT)), P.Opt, h.Opt, Color)
			}
		print("Optimal smoothing parameters:")
		print(paste(names(MyData)[i], "population:"))
		print(paste("Bandwidth: ", h.Opt*2+1, " observations"))
		print(paste("Polynomial degree: ", P.Opt))
		print(paste("Values of the AIC: ", round(PartOne[[i]]$AICMAT[h.Opt,P.Opt+1],4)))
		print(paste("Fitted DF: ", round(PartOne[[i]]$V2MAT[h.Opt,P.Opt+1],4)))
		M4B[[i]] <- c(FctMethod4_2ndPart(MyData[[i]]$Dxt, MyData[[i]]$Ext, MyData[[i]]$QxtRef, AgeRange, MyData[[i]]$AgeRef, MyData[[i]]$YearCom, MyData[[i]]$YearRef, P.Opt, h.Opt),list(AgeMethod=AgeRange,AgeCrit=PartOne[[1]]$AgeCrit,P.Opt=P.Opt,h.Opt=h.Opt))
	}
	if(Plot == T){
		.PlotMethod(M4B, MyData, AgeRange, Color)
	}	
	return(M4B)
}
