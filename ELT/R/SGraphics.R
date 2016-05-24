

## ------------------------------------------------------------------------ ##
##  Script R for  "Constructing Entity Specific Prospective Mortality Table ##
##                 Adjustment to a reference"                               ##
## ------------------------------------------------------------------------ ##
##  Script        SGraphics.R                                               ##
## ------------------------------------------------------------------------ ##
##  Description   Graphical Functions                                       ##
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
##  SurfacePlot function                                                    ##
## ------------------------------------------------------------------------ ##




SurfacePlot = function(xx, zexpr, mainexpr, axis, cc){
	XMat <- matrix(xx, nrow(xx), ncol(xx))
	if(ncol(XMat) < 10){ByNber = 1} else {ByNber = 5}
	jet.colors <- colorRampPalette(c("white", cc))
	color <- jet.colors(100)
	wireframe(XMat, zlab=list(zexpr, rot=90), ylab=list("Year"), xlab=list("Age"), main= mainexpr, aspect=c(2,1), drape = T, colorkey = F, screen = list(z = 20, x = -75), scales = list(arrows = F, y = list(at = c(seq(1,ncol(XMat), by = ByNber)), labels = c(seq(axis[3],axis[4],by=ByNber)), cex=.8), x = list( at = c(seq(1,nrow(XMat), by = 20)), labels = c(seq(axis[1],axis[2], by = 20)),cex=.8)), col.regions =  color) }
	
## ------------------------------------------------------------------------ ##
##  Graphics of the observed surfaces                                       ##
## ------------------------------------------------------------------------ ##

.PlotHistory = function(MyData, Color){
	YearBeg <- as.numeric(colnames(MyData[[1]]$Deaths)[1])
	YearEnd <- as.numeric(colnames(MyData[[1]]$Deaths)[ncol(MyData[[1]]$Deaths)])
	Path <- "Results/Graphics/Data"
	.CreateDirectory(paste("/",Path,sep=""))
	print(paste("Create the graphics of the observed statistics in .../",Path," ...", sep=""))
	for (i in 1 : (length(MyData)-1)){	
		png(filename=paste(Path,"/Indi",names(MyData)[i],".png", sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
		print(SurfacePlot(as.matrix(MyData[[i]]$Indi),expression(L[xt]), paste("Number of individuals,",names(MyData)[i],"pop."), c(0,120,YearBeg, YearEnd),Color))
		dev.off()
		png(filename=paste(Path,"/Expo",names(MyData)[i],".png", sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
		print(SurfacePlot(as.matrix(MyData[[i]]$Expo),expression(E[xt]), paste("Exposition,",names(MyData)[i],"pop."), c(0,120,YearBeg,YearEnd),Color))
		dev.off()
		png(filename=paste(Path,"/Deaths",names(MyData)[i],".png", sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
		print(SurfacePlot(as.matrix(MyData[[i]]$Deaths),expression(D[xt]), paste("Number of deaths,",names(MyData)[i],"pop."), c(0,120,YearBeg,YearEnd), Color))
		dev.off()
		png(filename=paste(Path,"/Proba",names(MyData)[i],".png", sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
		print(SurfacePlot(as.matrix(MyData[[i]]$Deaths/MyData[[i]]$Expo),expression(q[xt]), paste("Probability of death,",names(MyData)[i],"pop."), c(0,120,YearBeg,YearEnd), Color))
		dev.off()
		png(filename=paste(Path,"/LogProba",names(MyData)[i],".png", sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
		print(SurfacePlot(as.matrix(log(MyData[[i]]$Deaths/MyData[[i]]$Expo)),expression(paste("log ",q[xt])), paste("Probability of death,",names(MyData)[i],"pop. (Log)"), c(0,120,YearBeg,YearEnd), Color))
		dev.off()
		}
	}

## ------------------------------------------------------------------------ ##
##  Graphics of the extrapolated surfaces                                   ##
## ------------------------------------------------------------------------ ##

.PlotMethod = function(FitMethod, MyData, AgeMethod, Color){
	Path <- "Results/Graphics/Completion"
	.CreateDirectory(paste("/",Path,sep=""))
	print(paste("Create graphics of the fits before completion in .../",Path," ...", sep=""))
	for (i in 1 : length(FitMethod)){
		Age <- as.numeric(rownames(FitMethod[[i]]$QxtFitted))
		YearA <- min(MyData[[i]]$YearCom)
		YearB <- max(MyData[[i]]$YearRef)
		png(filename=paste(Path,"/",FitMethod[[i]]$NameMethod,"-BeforeCompletion-ProbaFitted",names(FitMethod)[i],".png", sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
		print(SurfacePlot(as.matrix(FitMethod[[i]]$QxtFitted),expression(widetilde(q)[xt]), paste(FitMethod[[i]]$NameMethod, "- Fitted prob. of death before completion,",names(MyData)[i],"pop."), c(min(Age), max(Age), YearA, YearB),Color))
		dev.off()
		png(filename=paste(Path,"/",FitMethod[[i]]$NameMethod,"-BeforeCompletion-LogProbaFitted",names(FitMethod)[i],".png", sep=""), width  = 2100, height = 2100, res=300, pointsize= 12)
		print(SurfacePlot(as.matrix(log(FitMethod[[i]]$QxtFitted)),expression(paste("log ",widetilde(q)[xt])), paste(FitMethod[[i]]$NameMethod, "- Fitted prob. of death before completion,",names(MyData)[i],"pop. (log)"), c(min(Age), max(Age), YearA, YearB), Color))
		dev.off()
		}
	}

## ------------------------------------------------------------------------ ##
##  Functions to plot the values of the criterion                           ##
## ------------------------------------------------------------------------ ##

.PlotCrit = function(crit, v2, LimCrit, LimV2) { 
	FigName <- paste("Values of the AIC for various polynomial degrees")
	VECA <- c(crit); VECB <- c(v2)
	xyplot(VECA ~ VECB|as.factor(FigName),xlab = "Fitted DF", ylab = "AIC", type = "n", par.settings = list(fontsize = list(text = 15, points = 12)), par.strip.text = list(cex = 1, lines = 1.1), strip = strip.custom(bg = "lightgrey"), scales = list(y = list(relation = "free", limits = list(LimCrit), at = list(round(c(seq(min(LimCrit),max(LimCrit), by = (max(LimCrit) - min(LimCrit)) / 10))),2)), x = list(relation = "free", limits = LimV2, tck = c(1, 1))), panel = function(...) {
		llines(VECA[1:20] ~ VECB[1:20],type = "p", pch = "0", col = 1)
		llines(VECA[21:40] ~ VECB[21:40],type = "p", pch = "1", col = 1)
		llines(VECA[41:60] ~ VECB[41:60],type = "p", pch = "2", col = 1)
		llines(VECA[61:80] ~ VECB[61:80],type = "p", pch = "3", col = 1) } ) }

.PlotCritChoice = function(crit, v2, LimCrit, LimV2, pp, hh, Color) { 
	FigName <- paste("Values of the AIC for various polynomial degrees")
	VECA <- c(crit); VECB <- c(v2)
	fig <- xyplot(VECA ~ VECB|as.factor(FigName),xlab = "Fitted DF", ylab= "AIC", type="n", par.settings = list(fontsize = list(text = 15, points = 12)), par.strip.text = list(cex = 1, lines = 1.1), strip = strip.custom(bg = "lightgrey"), scales = list(y = list(relation = "free", limits = list(LimCrit), at = list(round(c(seq(min(LimCrit),max(LimCrit),by = (max(LimCrit)-min(LimCrit))/10))),2)), x = list(relation = "free", limits = LimV2, tck = c(1,1))), panel = function(...) {
		llines(VECA[1:20]~VECB[1:20],type = "p", pch = "0", col = 1)
		llines(VECA[21:40]~VECB[21:40],type = "p", pch = "1", col = 1)
		llines(VECA[41:60]~VECB[41:60],type = "p", pch = "2", col = 1)
		llines(VECA[61:80]~VECB[61:80],type = "p", pch = "3", col = 1) } ) 
		print(fig)
trellis.focus("panel", 1, 1, highlight = FALSE)
panel.points(v2[hh,pp+1], crit[hh,pp+1], pch=16, col=Color, cex=3, lwd = 2)
panel.points(v2[hh,pp+1], crit[hh,pp+1], pch=as.character(pp), col="white", cex=1.5)		
trellis.unfocus()
		}

## ------------------------------------------------------------------------ ##
##  Plot Before/After                                                       ##
## ------------------------------------------------------------------------ ##	
.BeforeAfterCompletion = function(ModCompletion, OutputMethod, MyData, Age, Pop, Color, j){
	AgeObs <- (min(Age) - min(MyData$AgeRef) + 1) : nrow(MyData$Dxt)	
	LimMin <- min((log(MyData$Dxt/MyData$Ext)[AgeObs,as.character(j)])[log(MyData$Dxt/MyData$Ext)[AgeObs,as.character(j)] != -Inf], (log(OutputMethod$QxtFitted[,as.character(j)]))[log(OutputMethod$QxtFitted[,as.character(j)]) != -Inf],(log(ModCompletion$QxtFinal[,as.character(j)]))[log(ModCompletion$QxtFinal[,as.character(j)]) != -Inf] ,na.rm=T)
	LimMax <-  max((log(MyData$Dxt/MyData$Ext)[AgeObs,as.character(j)])[log(MyData$Dxt/MyData$Ext)[AgeObs,as.character(j)] != Inf], (log(OutputMethod$QxtFitted[,as.character(j)]))[log(OutputMethod$QxtFitted[,as.character(j)]) != Inf],(log(ModCompletion$QxtFinal[,as.character(j)]))[log(ModCompletion$QxtFinal[,as.character(j)]) != Inf],na.rm=T)
	print(plot(min(Age):130,c(log(MyData$Dxt/MyData$Ext)[AgeObs,as.character(j)],rep(NA,10)), ylim=c(LimMin,LimMax), type="p", pch=16, xlab="Age", ylab=expression(paste("log ",q[xt])), main=paste("Before/After Completion -", ModCompletion$NameMethod, "- Year",j,"-",Pop ,"pop.")))
	print(points(min(Age):(length(OutputMethod$QxtFitted[,as.character(j)])+min(Age)-1),log(OutputMethod$QxtFitted[,as.character(j)]),col=Color, type="l", lwd=2, lty=3))
	print(points(min(Age):130,log(ModCompletion$QxtFinal[,as.character(j)]),col=Color, type="l", lwd=2))
	legend(x=35,y=0,legend=c("Observations", "Before Completion","After Completion"), col=c("black",Color,Color), lty=c(NA,3,1),lwd=c(NA,2,2),pch=c(16,NA,NA), cex=.8, box.col="transparent")
	}

## ------------------------------------------------------------------------ ##
##  Plot of fits after completion for 2 populations                         ##
## ------------------------------------------------------------------------ ##
	
.FitPopsAfterCompletionLog = function(ModCompletion, MyData, Age, Color, j){
	AgeObs <- (min(Age) - min(MyData[[1]]$AgeRef) + 1) : nrow(MyData[[1]]$Dxt)
	LimMin <- vector(,2)
	for(i in 1:2){
		LimMin[i] <- min((log(MyData[[i]]$Dxt/MyData[[i]]$Ext)[AgeObs,as.character(j)])[log(MyData[[i]]$Dxt/MyData[[i]]$Ext)[AgeObs,as.character(j)] != -Inf],(log(ModCompletion[[i]]$QxtFinal[,as.character(j)]))[log(ModCompletion[[i]]$QxtFinal[,as.character(j)]) != -Inf] ,na.rm=T)
		}
	Lmin <- min(LimMin)
	print(plot(min(Age):130,c(log(MyData[[2]]$Dxt/MyData[[2]]$Ext)[AgeObs,as.character(j)],rep(NA,10)), ylim=c(Lmin,0.1), type="p", pch=16, xlab="Age", ylab=expression(paste("log ",q[xt])), main=paste("Fits after Completion -", ModCompletion[[1]]$NameMethod, "- Year",j)))
	print(points(min(Age):130,c(log(MyData[[1]]$Dxt/MyData[[1]]$Ext)[AgeObs,as.character(j)],rep(NA,10)), type="p", pch=1))
	print(points(min(Age):130,log(ModCompletion[[2]]$QxtFinal[,as.character(j)]),col=Color, type="l", lwd=2))
	print(points(min(Age):130,log(ModCompletion[[1]]$QxtFinal[,as.character(j)]),col=Color, type="l", lwd=2, lty=3))
	legend(x=35,y=0,legend=c(paste("Observations", names(MyData)[2]), paste("Observations", names(MyData)[1]), paste("Fit", names(MyData)[2]),paste("Fit", names(MyData)[1])), col=c(1,1,Color,Color), lty=c(NA,NA,1,3),lwd=c(NA,NA,2,2),pch=c(16,1,NA,NA), box.col="transparent",cex=.8)
	}

## ------------------------------------------------------------------------ ##
##  Plot of completion parameters                                           ##
## ------------------------------------------------------------------------ ##

.PlotParamCompletion = function(obs, Color) {
	k <- c(145,70)
	if(length(obs) == 2){
		Lim1 <- c(min(c(obs[[1]]$CompletionValues[,2], obs[[2]]$CompletionValues[,2])), max(c(obs[[1]]$CompletionValues[,2], obs[[2]]$CompletionValues[,2])))
		Lim2 <- c(min(c(obs[[1]]$CompletionValues[,3], obs[[2]]$CompletionValues[,3])), max(c(obs[[1]]$CompletionValues[,3], obs[[2]]$CompletionValues[,3])))
		NameObs <- names(obs)
		Pch <- c(1,16)
		Lty <- c(2,1)
		
	}
		if(length(obs) == 1){
		Lim1 <- c(min(obs[[1]]$CompletionValues[,2]), max(obs[[1]]$CompletionValues[,2]))
		Lim2 <- c(min(obs[[1]]$CompletionValues[,3]), max(obs[[1]]$CompletionValues[,3]))
		NameObs <- names(obs)
		Pch <- 1
		Lty <- 2
	}
	yy <- as.numeric(range(rownames(obs[[1]]$CompletionValues)))
	nn <- length(obs[[1]]$CompletionValues[,2])
	gr <- c(1,2)
	grp1 <- as.factor(rep(gr, each = nn))
	dat <- as.data.frame(cbind(c(obs[[1]]$CompletionValues[,2],obs[[1]]$CompletionValues[,3]), grp1))
	zzz <- with(dat, xyplot(c(obs[[1]]$CompletionValues[,2],obs[[1]]$CompletionValues[,3])  ~ rep(1:nn, length(gr)) | as.factor(grp1), 
	data = dat, layout = c(1,2), as.table = T, xlab = "Year", type="n",
	ylab = '', col = "black", scales = list(x = list(alternating = 1, 
	tck = c(1,0), at = c(seq(1,nn,5)), 
	labels = c(seq(yy[1],yy[2],5)) ), y = list(relation = "free",
	alternating = 3, tck = c(1,1), limits=list(c(Lim1[1],Lim1[2]), c(Lim2[1],Lim2[2])) )), par.settings = list(fontsize = list(text = 15,
	points = 7)), between=list(y = 0.6, x = 0.6), strip = strip.custom(which.given = 1,
	factor.levels = c(expression(paste( R^2," optimal")), expression(paste(c[t], " optimal"))), bg ="lightgrey"), par.strip.text = list(cex = 1,
	lines = 1.1))); print(zzz)
	trellis.focus("panel", 1, 1, highlight = FALSE)
		for (i in 1:length(obs)){
			panel.points(1:nrow(obs[[i]]$CompletionValues),obs[[i]]$CompletionValues[,2],type="b",lty=c(2,1)[i], pch=c(1,16)[i], col=Color)
		}
	trellis.unfocus()
	trellis.focus("panel", 1, 2, highlight = FALSE)
		for (i in 1:length(obs)){
			panel.points(1:nrow(obs[[i]]$CompletionValues),obs[[i]]$CompletionValues[,3],type="b",lty=c(2,1)[i], pch=c(1,16)[i], col=Color)
		}
	trellis.unfocus()
k1 <- draw.key(list(lines = list(pch = Pch, lty = Lty, type = "b", col = Color), cex = .9,
text = list(lab = NameObs), draw=TRUE))
pushViewport(viewport(x = unit(k[1], "mm"), y = unit(k[2], "mm"),
                 width = unit(1, "grobwidth", k1), 
                 height = unit(1, "grobheight", k1)))
                 grid.draw(k1);
	 }

## ------------------------------------------------------------------------ ##
##  Plot of the fits for the years in common original scale                 ##
## ------------------------------------------------------------------------ ##

.CompFittedyear = function(x, obs, yy ,k, ColVec, AgeVec, AgeRef, LimitY, pop, NameMethod) {
	yy <- as.character(yy)
	expr <- paste(NameMethod, pop,"- Comparison of the fits - year:", yy)
	xx <- matrix(,length(AgeVec),length(x))
	xx <- x[AgeVec+1-min(as.numeric(rownames(x))), yy]
	nn <- length(c(xx)); gr <- 1; grp1 <- as.factor(rep(gr, nn))
	dat <- as.data.frame(cbind(c(xx), grp1))
	zzz <- with(dat, xyplot(c(xx) ~ rep(AgeVec, length(gr)) | as.factor(grp1), 
	data = dat, layout = c(1,1), as.table = T, xlab = 'Age',
	ylab = expression(widetilde(q)[xt]), type= "n", scales = list(x = list(alternating = 1, tck = c(1,0),
	at = seq(min(AgeVec),max(AgeVec)+1,10), labels = seq(min(AgeVec),max(AgeVec)+1,10)), y = list(alternating = 3, 
	tck = c(1,1), limits= LimitY)), par.settings = list(fontsize = list(text = 15, points = 7)),
	between=list(y = 0.6, x = 0.6), par.strip.text = list(cex = 1, lines = 1.1),
	strip = strip.custom(which.given = 1, factor.levels = expr, bg = "lightgrey")))
	print(zzz)
	trellis.focus("panel", 1, 1, highlight = FALSE)
	panel.points(AgeVec, obs[AgeVec+1-min(AgeRef), yy], type="p", col=1, pch=16)
	panel.points(AgeVec, x[AgeVec+1-min(as.numeric(rownames(x))),yy], type="l", col=ColVec, lty = 1)
	trellis.unfocus()
	k1 <- draw.key(list(lines = list( type = c("p","l"), pch=16, col = c(1,ColVec), lty = 1), cex = .9, text = list(lab = c("Observations", NameMethod)), draw=TRUE))
	pushViewport(viewport(x = unit(k[1], "mm"), y = unit(k[2], "mm"), width = unit(1, "grobwidth", k1), height = unit(1, "grobheight", k1)))
	grid.draw(k1)
	}

.PlotFittedYear = function(OutputMethod, MyData, AgeCrit, Color, j, pop){
	LimMin <- min(((MyData$Dxt/MyData$Ext)[AgeCrit - min(MyData$AgeRef)+1,as.character(j)])[(MyData$Dxt/MyData$Ext)[AgeCrit - min(MyData$AgeRef)+1,as.character(j)] != -Inf], (OutputMethod$QxtFitted[AgeCrit-min(as.numeric(rownames(OutputMethod$QxtFitted)))+1,as.character(j)])[OutputMethod$QxtFitted[AgeCrit-min(as.numeric(rownames(OutputMethod$QxtFitted)))+1,as.character(j)] != -Inf],na.rm=T)
	LimMax <-  max(((MyData$Dxt/MyData$Ext)[AgeCrit-min(MyData$AgeRef)+1,as.character(j)])[(MyData$Dxt/MyData$Ext)[AgeCrit-min(MyData$AgeRef)+1,as.character(j)] != Inf], (OutputMethod$QxtFitted[AgeCrit-min(as.numeric(rownames(OutputMethod$QxtFitted)))+1,as.character(j)])[OutputMethod$QxtFitted[AgeCrit-min(as.numeric(rownames(OutputMethod$QxtFitted)))+1,as.character(j)] != Inf],na.rm=T)
	print(.CompFittedyear(OutputMethod$QxtFitted, MyData$Dxt/MyData$Ext, j , c(55,140), Color, AgeCrit, MyData$AgeRef, c(LimMin-0.01,LimMax+0.015), pop, OutputMethod$NameMethod))
	}

## ------------------------------------------------------------------------ ##
##  Plot of the fits for the years in common log scale                      ##
## ------------------------------------------------------------------------ ##

.CompFittedyearLog = function(x, obs, yy ,k, ColVec, AgeVec, AgeRef,  LimitY, pop, NameMethod) {
	yy <- as.character(yy)
	expr <- paste(NameMethod,pop,"- Comparison of the fits - year:", yy)
	xx <- x[AgeVec+1-min(as.numeric(rownames(x))), yy]
	nn <- length(c(xx)); gr <- 1; grp1 <- as.factor(rep(gr, nn))
	dat <- as.data.frame(cbind(c(xx), grp1))
	zzz <- with(dat, xyplot(c(xx) ~ rep(AgeVec, length(gr)) | as.factor(grp1), 
	data = dat, layout = c(1,1), as.table = T, xlab = 'Age',
	ylab = expression(paste("log ",widetilde(q)[xt])), type= "n", scales = list(x = list(alternating = 1, tck = c(1,0),
	at = seq(min(AgeVec),max(AgeVec)+1,10), labels = seq(min(AgeVec),max(AgeVec)+1,10)), y = list(alternating = 3, 
	tck = c(1,1), limits= LimitY)), par.settings = list(fontsize = list(text = 15, points = 7)),
	between=list(y = 0.6, x = 0.6), par.strip.text = list(cex = 1, lines = 1.1),
	strip = strip.custom(which.given = 1, factor.levels = expr, bg = "lightgrey")))
	print(zzz)
	trellis.focus("panel", 1, 1, highlight = FALSE)
	panel.points(AgeVec, log(obs[AgeVec+1-min(AgeRef), yy]), type="p", col=1, pch=16)
	panel.points(AgeVec, log(x[AgeVec+1-min(as.numeric(rownames(x))),yy]), type="l", col=ColVec, lty = 1)
	trellis.unfocus()
k1 <- draw.key(list(lines = list( type = c("p","l"), pch=16, col = c(1,ColVec), lty = 1), cex = .9, text = list(lab = c("Observations", NameMethod)), draw=TRUE))
	pushViewport(viewport(x = unit(k[1], "mm"), y = unit(k[2], "mm"), width = unit(1, "grobwidth", k1), height = unit(1, "grobheight", k1)))
	grid.draw(k1)
	}

.PlotFittedYearLog = function(OutputMethod, MyData, AgeCrit, Color, j, pop){
	LimMin <- min((log(MyData$Dxt/MyData$Ext)[AgeCrit - min(MyData$AgeRef)+1,as.character(j)])[log(MyData$Dxt/MyData$Ext)[AgeCrit - min(MyData$AgeRef)+1,as.character(j)] != -Inf], log((OutputMethod$QxtFitted[AgeCrit-min(as.numeric(rownames(OutputMethod$QxtFitted)))+1,as.character(j)]))[log(OutputMethod$QxtFitted[AgeCrit-min(as.numeric(rownames(OutputMethod$QxtFitted)))+1,as.character(j)]) != -Inf],na.rm=T)
	LimMax <-  max((log(MyData$Dxt/MyData$Ext)[AgeCrit-min(MyData$AgeRef)+1,as.character(j)])[log(MyData$Dxt/MyData$Ext)[AgeCrit-min(MyData$AgeRef)+1,as.character(j)] != Inf], log((OutputMethod$QxtFitted[AgeCrit-min(as.numeric(rownames(OutputMethod$QxtFitted)))+1,as.character(j)]))[log(OutputMethod$QxtFitted[AgeCrit-min(as.numeric(rownames(OutputMethod$QxtFitted)))+1,as.character(j)]) != Inf],na.rm=T)	
	print(.CompFittedyearLog(OutputMethod$QxtFitted, MyData$Dxt/MyData$Ext, j , c(55,140), Color, AgeCrit, MyData$AgeRef, c(LimMin-0.15,LimMax+0.15), pop,OutputMethod$NameMethod))
	}

## ------------------------------------------------------------------------ ##
##  Plot of the residuals                                                   ##
## ------------------------------------------------------------------------ ##

.PlotRes = function(Obs, xx, Color, yy, pop) {
	Names1 <- Obs$NameMethod
	Names2 <- names(Obs)[1:3]
	yy <- as.character(yy)
	x1 <- x2 <- x3 <- matrix(, length(xx), length(Names1))
	for (i in 1 : length(Names1)){
		x1[, i] <- Obs[[1]][, yy]
		x2[, i] <- Obs[[2]][, yy]
		x3[, i] <- Obs[[3]][, yy] }
	MAT <- matrix(NA,length(xx)*length(Names1),length(Names2))
	MAT[,1] <- c(x1)
	MAT[,2] <- c(x2)
	MAT[,3] <- c(x3)
	VEC <- c(MAT)
	grp1 <- rep(1:length(Names1), each = length(xx)*length(Names2))
	grp2 <- rep(1:length(Names2),times = length(Names1), each = length(xx))
	dat <- as.data.frame(cbind(VEC, grp1, grp2))
	zzz <- with(dat, xyplot(VEC ~ rep(xx,length(Names1)*length(Names2)) | as.factor(grp2) + 
as.factor(grp1), data = dat, as.table = T, between=list(y = 0.6, x = 0.6),
type = "n", main = paste(Names1,pop,"- Plot of the residuals - year:",yy),
xlab = 'Age', ylab = "", layout=c(length(Names1),length(Names2)),
par.settings = list(fontsize = list(text = 15, points = 7)),
par.strip.text = list(cex = .8, lines = 1.1),
scales = list(y = list(rot = 0, relation = "free", limits=c(list(c(range(x1[,1])),c(range(x2[,1])),c(range(x3[,1]))) )  ),
x = list(relation = "same", alternating = 1, tck = c(1,0), at = seq(min(xx),max(xx),20), labels = seq(min(xx),max(xx),20) )) ))

print(useOuterStrips(zzz, 
	strip.left = strip.custom(factor.levels = Names1, horizontal = F, bg = "lightgrey"), 
	strip = strip.custom(factor.levels = Names2, bg = "lightgrey"), 		strip.lines = 1.1, strip.left.lines = 1.1)) 
	for (i in 1:length(Names1)){
		trellis.focus("panel", 1, i, highlight = FALSE)
		panel.abline(h = 0, lty = 2, col = 1, cex=.5)
		panel.points(xx, x1[,i], type="p", pch=16, col=Color, cex=.5)
		fit.res <- locfit(x1[,i]  ~ xx, family = "gaussian", alpha = 0.15, link = "identity", kern="gauss")
		panel.points(xx, fitted(fit.res), type="l", col=Color, cex=.5)
		trellis.unfocus()
		trellis.focus("panel", 2, i, highlight = FALSE)
		panel.abline(h = 2, lty = 2, col = 1, cex=.5)
		panel.abline(h = -2, lty = 2, col = 1, cex=.5)
		panel.abline(h = 3, lty = 3, col = 1, cex=.5)
		panel.abline(h = -3, lty = 3, col = 1, cex=.5)
		panel.points(xx, x2[,i], type="p", pch=16, col=Color, cex=.5)
		fit.res2 <- locfit(x2[,i]  ~ xx, family = "gaussian", alpha = 45/95, link = "identity", kern="trwt")
		panel.points(xx, fitted(fit.res2), type="l", col=Color, cex=.5)
		trellis.unfocus()
		trellis.focus("panel", 3, i, highlight = FALSE)
		panel.abline(h = 0, lty = 2, col = 1, cex=.5)
        panel.points(xx, x3[,i], type = "b",col = Color, pch=16,cex=.5)
		trellis.unfocus()
		}
	}

## ------------------------------------------------------------------------ ##
##  Pointwise confidence intervals of the fitted deaths                     ##
## ------------------------------------------------------------------------ ##

.PlotDIntConf = function(d, obs, AgeVec, yy , ColVec, pop) {
	x <- d$DxtFitted; AgeRef <- obs$AgeRef; k <- c(55,140)
	yy <- as.character(yy)
	expr <- paste(d$NameMethod,pop,"- Comparison fitted deaths - year:", yy)
	xx <- x[AgeVec+1-min(as.numeric(rownames(x))), yy]
	nn <- length(c(xx)); gr <- 1; grp1 <- as.factor(rep(gr, nn))
	dat <- as.data.frame(cbind(c(xx), grp1))
	zzz <- with(dat, xyplot(c(xx) ~ rep(AgeVec, length(gr)) | as.factor(grp1), 
	data = dat, layout = c(1,1), as.table = T, xlab = 'Age',
	ylab = expression(widetilde(D)[xt]), type= "n", scales = list(x = list(alternating = 1, tck = c(1,0),
	at = seq(min(AgeVec),max(AgeVec)+1,10), labels = seq(min(AgeVec),max(AgeVec)+1,10)), y = list(alternating = 3, 
	tck = c(1,1), limits= c(-1,max(unlist(d$DIntUp))+1))), par.settings = list(fontsize = list(text = 15, points = 7)),
	between=list(y = 0.6, x = 0.6), par.strip.text = list(cex = 1, lines = 1.1),
	strip = strip.custom(which.given = 1, factor.levels = expr, bg = "lightgrey")))
	print(zzz)
	trellis.focus("panel", 1, 1, highlight = FALSE)
	panel.points(AgeVec, obs$Dxt[AgeVec+1-min(AgeRef), yy], type="p", col=1, pch=16)
		panel.points(AgeVec, x[AgeVec+1-min(as.numeric(rownames(x))),yy], type="l", col=ColVec, lty = 1)
		panel.points(AgeVec, d$DIntUp[AgeVec+1-min(as.numeric(rownames(d$DIntUp))),yy], type="l", col=ColVec, lty = 3)
		panel.points(AgeVec, d$DIntLow[AgeVec+1-min(as.numeric(rownames(d$DIntLow))),yy], type="l", col=ColVec, lty = 3)
	trellis.unfocus()
k1 <- draw.key(list(lines = list( type = c("p","l"), pch=16, col = c(1,ColVec), lty = 1), cex = .9, text = list(lab = c("Observations", d$NameMethod)), draw=TRUE))
	pushViewport(viewport(x = unit(k[1], "mm"), y = unit(k[2], "mm"), width = unit(1, "grobwidth", k1), height = unit(1, "grobheight", k1)))
	grid.draw(k1)
	}

## ------------------------------------------------------------------------ ##
##  Plot of the evolution of the periodic life expectancy                   ##
## ------------------------------------------------------------------------ ##

.PlotPerExp = function(x1, x2, x3, aa, aaa,k, col.vec, year, lly, pop, leg.vec) {
	expr <- paste("Comparison of the trends -",leg.vec,"-",pop ,"pop. - age", aa)
	xx <- x2[aa+1-min(aaa),]
	nn <- length(c(xx)); Level <- 1; Levelp1 <- as.factor(rep(Level, nn))
	dat <- as.data.frame(cbind(c(xx), Levelp1))
	zzz <- with(dat, xyplot(c(xx) ~ rep(1:ncol(x2), length(Level)) | as.factor(Levelp1), 
	data = dat, layout = c(1,1), as.table = T, xlab = 'Year', 
	ylab = "Periodic life expectancies", type= "n", scales = list(x = list(alternating = 1, tck = c(1,0), cex=0.8,
	at = seq(1,ncol(x2),5), labels = seq(min(year),max(year),5)), y = list(alternating = 3, 
	tck = c(1,1), limits= lly)), par.settings = list(fontsize = list(text = 15, points = 7)),
	between=list(y = 0.6, x = 0.6), par.strip.text = list(cex = 1, lines = 1.1),
	strip = strip.custom(which.given = 1, factor.levels = expr, bg = "lightgrey")))
	print(zzz)
	trellis.focus("panel", 1, 1, highlight = FALSE)
	panel.points(1:ncol(x1), x1[aa+1-min(aaa),], type="b", col=1, lty = 2, pch=16)
	panel.points(ncol(x1):ncol(x2), x2[aa+1-min(aaa),ncol(x1):ncol(x2)], type="l", col=col.vec, lty = 1)
	trellis.unfocus()
k1 <- draw.key(list(lines = list( type = c("b","l"), pch=16, col = c(1,col.vec), lty = c(2,1)), cex = .9, text = list(lab = c("Observations", leg.vec)), draw=TRUE))

pushViewport(viewport(x = unit(k[1], "mm"), y = unit(k[2], "mm"), width = unit(1, "grobwidth", k1), height = unit(1, "grobheight", k1)))

grid.draw(k1) }

## ------------------------------------------------------------------------ ##
##  Graphic of the fits of the simulated data                               ##
## ------------------------------------------------------------------------ ##

.SimPlot = function(q1, q2, MyData, x1, yy, pop, nsim, Color){
	if(length(MyData) == 3){
		LayOut <- c(2,1)
		Lim <- c(list(c(min(unlist(q2[[1]])),.15)), list(c(min(unlist(q2[[2]])),.15)))
		}
	if(length(MyData) == 2){
		LayOut <- c(1,1)
		Lim <- c(list(c(min(unlist(q2[[1]])),.15)))
		}	
	aa <- x1+1-min(x1)
	qq <- (MyData[[1]]$Dxt/MyData[[1]]$Ext)[x1+1-min(MyData[[1]]$AgeRef),as.character(yy)]
	MAT.A <- c(rep(q1[[1]]$QxtFinal[aa,as.character(yy)], length(pop)))
	VEC.A <- c(MAT.A)
	grp1 <- rep(pop, each = length(aa))
	dat <- as.data.frame(cbind(VEC.A, grp1))
	zzz <- with(dat, xyplot(VEC.A ~ rep(x1,length(pop)) | as.factor(grp1), data = dat, as.table = T, between=list(y = 0.6, x = 0.6),
type = "n", xlab = 'Age', ylab = "", layout=LayOut, par.settings = list(fontsize = list(text = 15, points = 7)), par.strip.text = list(cex = 1, lines = 1.1),
	strip = strip.custom(which.given = 1, factor.levels = pop, bg = "lightgrey"), scales = list(y = list(relation = "free", rot = 0, limits = Lim), x = list(relation = "same", draw=T, alternating = 1, tck = c(1,0), limits =  c(min(x1)-2, max(x1+2))))))
	print(zzz)
	for(i in 1 : (length(MyData)-1)){
    trellis.focus("panel", i, 1, highlight = FALSE)
       for (k in 1:nsim) {
		panel.points(x1, q2[[i]][[k]][aa,as.character(yy)],type="l",col=Color) }
		panel.points(x1, (MyData[[i]]$Dxt/MyData[[i]]$Ext)[x1+1-min(MyData[[i]]$AgeRef),as.character(yy)],type="p",col=1, pch=16)
		panel.points(x1, q1[[i]]$QxtFinal[aa,as.character(yy)],type="l",col=1, lwd=2)
    trellis.unfocus()
		}
	}

## ------------------------------------------------------------------------ ##
##  Plot quantile of the simulated life expectancies                        ##
## ------------------------------------------------------------------------ ##

.PlotExpQtle = function(x, yy , age.vec, expr, cc) {
	yyy <- seq(min(yy),max(yy),by=10)
	yyy <- as.character(yyy)
	av <- seq(min(age.vec),max(age.vec),by=10)
	ly <- length(yyy)*2 + 1
	maxx <- minx <- vector(,length(yyy))
	for(i in 1:length(yyy)){
		maxx[i] <- x[[which(yy==yyy[i])]][min(av)-min(age.vec)+1,2]
		minx[i] <- x[[which(yy==yyy[i])]][max(av)-min(age.vec)+1,2]
	}
	maxxx <- ceiling(max(maxx)); minxx <- floor(min(minx))
	lx <- seq(minxx,maxxx,length.out=ly)
	nn <- ly; gr <- 1; grp1 <- as.factor(rep(gr, nn))
	dat <- as.data.frame(cbind(c(ly), grp1))
	zzz <- with(dat, xyplot(c(1:ly) ~ lx | as.factor(grp1), 
	data = dat, layout = c(1,1), as.table = T, xlab = 'Periodic life expectancy',
	ylab = "", type= "n", scales = list(x = list(alternating = 1, tck = c(1,0),
	at = seq(minxx,maxxx,10), labels = seq(minxx,maxxx,10)), y = list(alternating = 3, 
	tck = c(1,1),at = seq(2.5,ly,by=2), labels = yyy)), par.settings = list(fontsize = list(text = 15, points = 7)),
	between=list(y = 0.6, x = 0.6), par.strip.text = list(cex = 1, lines = 1.1),
	strip = strip.custom(which.given = 1, factor.levels = expr, bg = "lightgrey")))
	print(zzz)
	trellis.focus("panel", 1, 1, highlight = FALSE)
	rr <- seq(2,length(yyy)*2,by=2)
	for(i in 1:length(av)){
		aa <- vector(,length(yyy))
		for(j in 1:length(yyy)){
			tt <- x[[which(yy==yyy[j])]][av[i]-min(age.vec)+1,]
			panel.rect(tt[1], rr[j], tt[3], rr[j]+1,col=cc, border="black") 
			for(k in 1:3){
				panel.points(rep(tt[k], 2),c(rr[j],rr[j]+1), lwd=2, col=1, type="l")
				}
			aa[j] <- tt[2]
			}
		panel.text(mean(aa), 1.25,labels=paste("Age ",av[i]))
		}
	trellis.unfocus() 
	}

## ------------------------------------------------------------------------ ##
##  Plot of the relative dispersion                                         ##
## ------------------------------------------------------------------------ ##

.PlotRelDisp = function(xx, yy , age.vec, expr, cc, k) {
	yyy <- seq(min(yy),max(yy),by=10)
	yyy <- as.character(yyy)
	av <- seq(min(age.vec),max(age.vec),by=10)
	ly <- length(yyy)
	xxx <- matrix(,length(av),length(yyy))
	for(j in 1:length(yyy)){
		xxx[,j] <- xx[av-min(age.vec)+1, which(yy==yyy[j])] * 100
	}
	minxx <- floor((min(unlist(xxx), na.rm=T)*100))/100
	maxx <- ceiling((max(unlist(xxx), na.rm=T)*100))/100
	lx <- seq(minxx,maxx,length.out=ly)
	nn <- ly; gr <- 1; grp1 <- as.factor(rep(gr, nn))
	dat <- as.data.frame(cbind(c(ly), grp1))
	zzz <- with(dat, xyplot(c(1:ly) ~ lx | as.factor(grp1), 
	data = dat, layout = c(1,1), as.table = T, xlab = 'Relative dispersion (%)',
	ylab = "", type= "n", scales = list(x = list(alternating = 1, tck = c(1,0),
	at = seq(minxx,maxx,length.out=10), labels = round(seq(minxx,maxx,length.out=10),2)), y = list(alternating = 3, 
	tck = c(1,1),at = seq(1,ly,by=1), labels = yyy, limits=c(0,7))), par.settings = list(fontsize = list(text = 15, points = 7)),
	between=list(y = 0.6, x = 0.6), par.strip.text = list(cex = 1, lines = 1.1),
	strip = strip.custom(which.given = 1, factor.levels = expr, bg = "lightgrey")))
	print(zzz)
	trellis.focus("panel", 1, 1, highlight = FALSE)
	ppp <- c(0:6,8)
		for(j in 1:length(av)){
			panel.points(xxx[j,],1:ly , lwd=2, col=cc, type="b", pch=ppp[j], lty =2, cex=1.5)
			}
		trellis.unfocus()
		k1 <- draw.key(list(lines = list( type = "b", pch=ppp, col = cc, lty = 2, cex=1.5), cex = 1.1, text = list(lab = paste("Age", av)), draw=TRUE))

pushViewport(viewport(x = unit(k[1], "mm"), y = unit(k[2], "mm"), width = unit(1, "grobwidth", k1), height = unit(1, "grobheight", k1)))

grid.draw(k1) }

## ------------------------------------------------------------------------ ##
##  Comparison of the fits - original scale                                 ##
## ------------------------------------------------------------------------ ##

.ComparisonFitsMethods = function(ListOutputs, j, MyData, AgeCrit, Pop, ColorComp, LtyComp){
	NamesOfMethods <- ColVec <- LtyVec <- vector(,length(ListOutputs))
	for(i in 1:length(ListOutputs)){
		NamesOfMethods[i] <- ListOutputs[[i]][[Pop]]$NameMethod
		if(NamesOfMethods[i] == "Method1"){
			ColVec[i] <- ColorComp[1]
			LtyVec[i] <- LtyComp[1]
			}
		if(NamesOfMethods[i] == "Method2"){
			ColVec[i] <- ColorComp[2]
			LtyVec[i] <- LtyComp[2]
			}
		if(NamesOfMethods[i] == "Method3"){
			ColVec[i] <- ColorComp[3]
			LtyVec[i] <- LtyComp[3]
			}
		if(NamesOfMethods[i] == "Method4"){
			ColVec[i] <- ColorComp[4]
			LtyVec[i] <- LtyComp[4]
			}
	}
	QxtFittedList <- vector("list",length(ListOutputs))
	names(QxtFittedList) <- NamesOfMethods
	for(i in 1:length(ListOutputs)){
		QxtFittedList[[i]] <- ListOutputs[[i]][[Pop]]$QxtFitted
		}
		.PlotFitsMethods(QxtFittedList, MyData, AgeCrit, ColVec, LtyVec, j, NamesOfMethods, Pop)
	}

.PlotFitsMethods = function(QxtFittedList, MyData, AgeCrit, ColVec, LtyVec, j, NamesOfMethods, pop){
	QxtFittedMethods <- vector("list",length(QxtFittedList))
	for(i in 1:length(QxtFittedList)){
		QxtFittedMethods[[i]] <- QxtFittedList[[i]][AgeCrit - min(as.numeric(rownames(QxtFittedList[[i]]))) + 1, as.character(j)]
	}
	LimMin <- min(((MyData$Dxt/MyData$Ext)[AgeCrit - min(MyData$AgeRef)+1,as.character(j)])[(MyData$Dxt/MyData$Ext)[AgeCrit - min(MyData$AgeRef)+1,as.character(j)] != -Inf], (unlist(QxtFittedMethods))[(unlist(QxtFittedMethods)) != -Inf],na.rm=T)
	LimMax <-  max(((MyData$Dxt/MyData$Ext)[AgeCrit - min(MyData$AgeRef)+1,as.character(j)])[(MyData$Dxt/MyData$Ext)[AgeCrit - min(MyData$AgeRef)+1,as.character(j)] != -Inf], (unlist(QxtFittedMethods))[(unlist(QxtFittedMethods))!= -Inf],na.rm=T)	
	print(.CompMethodsByYears(QxtFittedList, MyData$Dxt/MyData$Ext, j , c(55,140), ColVec, LtyVec, NamesOfMethods, AgeCrit, MyData$AgeRef, c(LimMin-0.01,LimMax+0.02), pop))
	}

.CompMethodsByYears = function(x, obs, yy ,k, col.vec, lty.vec, leg.vec, AgeCrit, ageref, lly, pop) {
	yy <- as.character(yy)
	expr <- paste("Comparison of the fits for the year: ", yy)
	xx <- log(obs[AgeCrit+1-min(ageref), yy])
	nn <- length(c(xx)); gr <- 1; grp1 <- as.factor(rep(gr, nn))
	dat <- as.data.frame(cbind(c(xx), grp1))
	zzz <- with(dat, xyplot(c(xx) ~ rep(AgeCrit, length(gr)) | as.factor(grp1), data = dat, layout = c(1,1), as.table = T, xlab = 'Age', ylab = expression(widetilde(q)[xt]), type= "n", scales = list(x = list(alternating = 1, tck = c(1,0), at = seq(min(AgeCrit),max(AgeCrit)+1,10), labels = seq(min(AgeCrit),max(AgeCrit)+1,10)), y = list(alternating = 3, tck = c(1,1), limits= lly)), par.settings = list(fontsize = list(text = 15, points = 7)), between=list(y = 0.6, x = 0.6), par.strip.text = list(cex = 1, lines = 1.1), strip = strip.custom(which.given = 1, factor.levels = expr, bg = "lightgrey")))
	print(zzz)
	trellis.focus("panel", 1, 1, highlight = FALSE)
	panel.points(AgeCrit, obs[AgeCrit+1-min(ageref), yy], type="p", col=1, pch=16)
	for( i in 1:length(x)) {
		panel.points(AgeCrit, x[[i]][AgeCrit+1-min(as.numeric(rownames(x[[i]]))),yy], type="l", col=col.vec[i], lty=lty.vec[i])
		}
	trellis.unfocus()
	k1 <- draw.key(list(lines = list( type = c("p",rep("l",length(x))), pch=16, col = c(1,col.vec), lty = c(NA,lty.vec)), cex = .9, text = list(lab = c("Observations", leg.vec)), draw=TRUE))
	pushViewport(viewport(x = unit(k[1], "mm"), y = unit(k[2], "mm"), width = unit(1, "grobwidth", k1), height = unit(1, "grobheight", k1)))
	grid.draw(k1)
	}

## ------------------------------------------------------------------------ ##
##  Comparison of the fits - log scale                                      ##
## ------------------------------------------------------------------------ ##

.ComparisonFitsMethodsLog = function(ListOutputs, j, MyData, AgeCrit, Pop, ColorComp, LtyComp){
	NamesOfMethods <- ColVec <- LtyVec <- vector(,length(ListOutputs))
	AgeList <- vector("list",length(ListOutputs))
	for(i in 1:length(ListOutputs)){
		NamesOfMethods[i] <- ListOutputs[[i]][[Pop]]$NameMethod
		if(NamesOfMethods[i] == "Method1"){
			ColVec[i] <- ColorComp[1]
			LtyVec[i] <- LtyComp[1]
			}
		if(NamesOfMethods[i] == "Method2"){
			ColVec[i] <- ColorComp[2]
			LtyVec[i] <- LtyComp[2]
			}
		if(NamesOfMethods[i] == "Method3"){
			ColVec[i] <- ColorComp[3]
			LtyVec[i] <- LtyComp[3]
			}
		if(NamesOfMethods[i] == "Method4"){
			ColVec[i] <- ColorComp[4]
			LtyVec[i] <- LtyComp[4]
			}
	}
	QxtFittedList <- vector("list",length(ListOutputs))
	names(QxtFittedList) <- NamesOfMethods
	for(i in 1:length(ListOutputs)){
		QxtFittedList[[i]] <- ListOutputs[[i]][[Pop]]$QxtFitted
		}
		.PlotFitsMethodsLog(QxtFittedList, MyData, AgeCrit, ColVec, LtyVec, j, NamesOfMethods, Pop)
	}

.PlotFitsMethodsLog = function(QxtFittedList, MyData, AgeCrit, ColVec, LtyVec, j, NamesOfMethods, pop){
	QxtFittedMethods <- vector("list",length(QxtFittedList))
	for(i in 1:length(QxtFittedList)){
		QxtFittedMethods[[i]] <- QxtFittedList[[i]][AgeCrit - min(as.numeric(rownames(QxtFittedList[[i]]))) + 1, as.character(j)]
	}
	LimMin <- min((log(MyData$Dxt/MyData$Ext)[AgeCrit - min(MyData$AgeRef)+1,as.character(j)])[log(MyData$Dxt/MyData$Ext)[AgeCrit - min(MyData$AgeRef)+1,as.character(j)] != -Inf], log(unlist(QxtFittedMethods))[log(unlist(QxtFittedMethods)) != -Inf],na.rm=T)
	LimMax <-  max((log(MyData$Dxt/MyData$Ext)[AgeCrit - min(MyData$AgeRef)+1,as.character(j)])[log(MyData$Dxt/MyData$Ext)[AgeCrit - min(MyData$AgeRef)+1,as.character(j)] != -Inf], log(unlist(QxtFittedMethods))[log(unlist(QxtFittedMethods))!= -Inf],na.rm=T)	
	print(.CompMethodsByYearsLog(QxtFittedList, MyData$Dxt/MyData$Ext, j , c(55,140), ColVec, LtyVec, NamesOfMethods, AgeCrit, MyData$AgeRef, c(LimMin-0.15,LimMax+0.15), pop))
	}

.CompMethodsByYearsLog = function(x, obs, yy ,k, col.vec, lty.vec, leg.vec, AgeCrit, ageref, lly, pop) {
	yy <- as.character(yy)
	expr <- paste(pop,"- Comparison of the fits for the year:", yy)
	xx <- log(obs[AgeCrit+1-min(ageref), yy])
	nn <- length(c(xx)); gr <- 1; grp1 <- as.factor(rep(gr, nn))
	dat <- as.data.frame(cbind(c(xx), grp1))
	zzz <- with(dat, xyplot(c(xx) ~ rep(AgeCrit, length(gr)) | as.factor(grp1), data = dat, layout = c(1,1), as.table = T, xlab = 'Age', ylab = expression(paste("log ", widetilde(q)[xt])), type= "n", scales = list(x = list(alternating = 1, tck = c(1,0), at = seq(min(AgeCrit),max(AgeCrit)+1,10), labels = seq(min(AgeCrit),max(AgeCrit)+1,10)), y = list(alternating = 3, tck = c(1,1), limits= lly)), par.settings = list(fontsize = list(text = 15, points = 7)), between=list(y = 0.6, x = 0.6), par.strip.text = list(cex = 1, lines = 1.1), strip = strip.custom(which.given = 1, factor.levels = expr, bg ="lightgrey")))
	print(zzz)
	trellis.focus("panel", 1, 1, highlight = FALSE)
	panel.points(AgeCrit, log(obs[AgeCrit+1-min(ageref), yy]), type="p", col=1, pch=16)
	for( i in 1:length(x)) {
		panel.points(AgeCrit, log(x[[i]][AgeCrit+1-min(as.numeric(rownames(x[[i]]))),yy]), type="l", col=col.vec[i], lty = lty.vec[i])
		}
	trellis.unfocus()
	k1 <- draw.key(list(lines = list( type = c("p",rep("l",length(x))), pch=16, col = c(1,col.vec), lty = c(NA, lty.vec)), cex = .9, text = list(lab = c("Observations", leg.vec)), draw=TRUE))
	pushViewport(viewport(x = unit(k[1], "mm"), y = unit(k[2], "mm"), width = unit(1, "grobwidth", k1), height = unit(1, "grobheight", k1)))
	grid.draw(k1)
	}

## ------------------------------------------------------------------------ ##
##  Comparison of the residuals                                             ##
## ------------------------------------------------------------------------ ##

.ComparisonResidualsMethods = function(ListValidationLevel1, j, AgeCrit, Pop, ColorComp ){
	NamesOfMethods <- ColVec <- vector(,length(ListValidationLevel1))
	for(i in 1:length(ListValidationLevel1)){
		NamesOfMethods[i] <- ListValidationLevel1[[i]][[1]][[1]]$NameMethod
		if(NamesOfMethods[i] == "Method1"){
			ColVec[i] <- ColorComp[1]
			}
		if(NamesOfMethods[i] == "Method2"){
			ColVec[i] <- ColorComp[2]
			}
		if(NamesOfMethods[i] == "Method3"){
			ColVec[i] <- ColorComp[3]
			}
		if(NamesOfMethods[i] == "Method4"){
			ColVec[i] <- ColorComp[4]
			}
	}
	RespResList <- PearResList <- DevResList <- vector("list",length(ListValidationLevel1))
	names(RespResList) <- names(PearResList) <- names(DevResList) <- NamesOfMethods
	for(i in 1:length(ListValidationLevel1)){
		RespResList[[i]] <- ListValidationLevel1[[i]][[2]][[Pop]]$"Response Residuals"
		PearResList[[i]] <- ListValidationLevel1[[i]][[2]][[Pop]]$"Pearson Residuals"				
		DevResList[[i]] <- ListValidationLevel1[[i]][[2]][[Pop]]$"Deviance Residuals"
		}
	print(.PlotResidualsMethods(RespResList,PearResList,DevResList, j, AgeCrit, NamesOfMethods, names(ListValidationLevel1[[1]][[2]][[Pop]])[1:3], ColVec, Pop))
	}


.PlotResidualsMethods = function(obs.1,obs.2,obs.3, yy, xx, names1, names2, colvec, pop) {
	yy <- as.character(yy)
	Lim <- vector("list", length(names1))	
	for(i in 1:length(names1)){
		Lim[[i]] <- list(c(range(obs.1[[i]][,yy])),c(range(obs.2[[i]][,yy])),c(range(obs.3[[i]][,yy])))
		}
	if(length(names1)==1){ LLim <- c(Lim[[1]]) }
	if(length(names1)==2){ LLim <- c(Lim[[1]], Lim[[2]]) }
	if(length(names1)==3){ LLim <- c(Lim[[1]], Lim[[2]], Lim[[3]]) }
	if(length(names1)==4){ LLim <- c(Lim[[1]], Lim[[2]], Lim[[3]], Lim[[4]]) }
	MAT <- matrix(NA,length(xx)*length(names1),length(names2))
	MAT[,1] <- rep(obs.1[[1]][,yy], length(names1))
	MAT[,2] <- rep(obs.1[[2]][,yy], length(names1))
	MAT[,3] <- rep(obs.1[[3]][,yy], length(names1))
	VEC <- c(MAT)
	grp1 <- rep(1:length(names1), each = length(xx)*length(names2))
	grp2 <- rep(1:length(names2),times = length(names1), each = length(xx))
	dat <- as.data.frame(cbind(VEC, grp1, grp2))
	zzz <- with(dat, xyplot(VEC ~ rep(xx,length(names1)*length(names2)) | as.factor(grp2) + 
as.factor(grp1), data = dat, as.table = T, between=list(y = 0.6, x = 0.6),
type = "n", main = paste(pop,"- Plot of the residuals - year:",yy),
xlab = 'Age', ylab = "", layout=c(length(names1),length(names2)),
par.settings = list(fontsize = list(text = 15, points = 7)),
par.strip.text = list(cex = .8, lines = 1.1),
scales = list(y = list(rot = 0, relation = "free", limits=c(LLim)  ),
x = list(relation = "same", alternating = 1, tck = c(1,0), at = seq(min(xx),max(xx),20), labels = seq(min(xx),max(xx),20) )) ))
	print(useOuterStrips(zzz, 
	strip.left = strip.custom(factor.levels = names1, horizontal = F, bg ="lightgrey"), 
	strip = strip.custom(factor.levels = names2, bg ="lightgrey"), 		strip.lines = 1.1, strip.left.lines = 1.1)) 
	for (i in 1:length(names1)){
		trellis.focus("panel", 1, i, highlight = FALSE)
		panel.abline(h = 0, lty = 2, col = 1, cex=.5)
		panel.points(xx, obs.1[[i]][,yy], type="p", pch=16, col=colvec[i], cex=.5)
		fit.res <- locfit(obs.1[[i]][,yy]  ~ xx, family = "gaussian", alpha = 0.15, link = "identity", kern="gauss")
		panel.points(xx, fitted(fit.res), type="l", col=colvec[i], cex=.5)
		trellis.unfocus()
		trellis.focus("panel", 2, i, highlight = FALSE)
		panel.abline(h = 2, lty = 2, col = 1, cex=.5)
		panel.abline(h = -2, lty = 2, col = 1, cex=.5)
		panel.abline(h = 3, lty = 3, col = 1, cex=.5)
		panel.abline(h = -3, lty = 3, col = 1, cex=.5)
		panel.points(xx, obs.2[[i]][,yy], type="p", pch=16, col=colvec[i], cex=.5)
		fit.res2 <- locfit(obs.2[[i]][,yy]  ~ xx, family = "gaussian", alpha = 45/95, link = "identity", kern="trwt")
		panel.points(xx, fitted(fit.res2), type="l", col=colvec[i], cex=.5)
		trellis.unfocus()
		trellis.focus("panel", 3, i, highlight = FALSE)
		panel.abline(h = 0, lty = 2, col = 1, cex=.5)
       	panel.points(xx, obs.3[[i]][,yy], type = "b",col = colvec[i], pch=16,cex=.5)
		trellis.unfocus()
		}
	}
	
## ------------------------------------------------------------------------ ##
##  Comparison of the fitted deaths                                         ##
## ------------------------------------------------------------------------ ##

.ComparisonFittedDeathsMethods = function(ListValidationLevel1, j, MyData, AgeCrit, Pop, ColorComp, LtyComp){
	NamesOfMethods <- ColVec <- LtyVec <- vector(,length(ListValidationLevel1))
	for(i in 1:length(ListValidationLevel1)){
		NamesOfMethods[i] <- ListValidationLevel1[[i]][[1]][[1]]$NameMethod
		if(NamesOfMethods[i] == "Method1"){
			ColVec[i] <- ColorComp[1]
			LtyVec[i] <- LtyComp[1]
			}
		if(NamesOfMethods[i] == "Method2"){
			ColVec[i] <- ColorComp[2]
			LtyVec[i] <- LtyComp[2]
			}
		if(NamesOfMethods[i] == "Method3"){
			ColVec[i] <- ColorComp[3]
			LtyVec[i] <- LtyComp[3]
			}
		if(NamesOfMethods[i] == "Method4"){
			ColVec[i] <- ColorComp[4]
			LtyVec[i] <- LtyComp[4]
			}
	}
	DxtFittedList <- DxtUpList <- DxtLowList <- vector("list",length(ListValidationLevel1))
	names(DxtFittedList) <- names(DxtUpList) <- names(DxtLowList) <- NamesOfMethods
	for(i in 1:length(ListValidationLevel1)){
		DxtFittedList[[i]] <- ListValidationLevel1[[i]][[3]][[Pop]]$DxtFitted
		DxtUpList[[i]] <- ListValidationLevel1[[i]][[3]][[Pop]]$DIntUp				
		DxtLowList[[i]] <- ListValidationLevel1[[i]][[3]][[Pop]]$DIntLow
		}
	print(.PlotFitteddeathsMethods(DxtFittedList, MyData$Dxt, DxtUpList, DxtLowList, j , c(55,140), ColVec, LtyVec, NamesOfMethods, AgeCrit, MyData$AgeRef, Pop))
	}

.PlotFitteddeathsMethods = function(x, obs, dup, dlow, yy ,k, col.vec, lty.vec, leg.vec, age.vec, ageref, Pop) {
	yy <- as.character(yy)
	expr <- paste(Pop,",comparison of the fitted deaths for the year:", yy)
	Lim <-  c(-2, max(unlist(dup))+2)
	nn <- length(c(x[[1]])); gr <- 1; grp1 <- as.factor(rep(gr, nn))
	dat <- as.data.frame(cbind(c(x[[1]]), grp1))
	zzz <- with(dat, xyplot(c(x[[1]]) ~ rep(age.vec, length(gr)) | as.factor(grp1), 
	data = dat, layout = c(1,1), as.table = T, xlab = 'Age',
	ylab = expression(widetilde(D)[xt]), type= "n", scales = list(x = list(alternating = 1, tck = c(1,0),
	at = seq(min(age.vec),max(age.vec)+1,10), labels = seq(min(age.vec),max(age.vec)+1,10)), y = list(alternating = 3, 
	tck = c(1,1), limits= Lim)), par.settings = list(fontsize = list(text = 15, points = 7)),
	between=list(y = 0.6, x = 0.6), par.strip.text = list(cex = 1, lines = 1.1),
	strip = strip.custom(which.given = 1, factor.levels = expr, bg ="lightgrey")))
	print(zzz)
	trellis.focus("panel", 1, 1, highlight = FALSE)
	panel.points(age.vec, obs[age.vec+1-min(ageref), yy], type="p", col=1, pch=16)
	if(col.vec[1] != col.vec[2]){
	for( i in 1:length(x)) { 
		panel.points(age.vec, x[[i]][age.vec+1-min(age.vec),yy], type="l", col=col.vec[i], lty = 1)
		panel.points(age.vec, dup[[i]][age.vec+1-min(age.vec),yy], type="l", col=col.vec[i], lty = 2)
		panel.points(age.vec, dlow[[i]][age.vec+1-min(age.vec),yy], type="l", col=col.vec[i], lty = 2)
		}
	trellis.unfocus()
	k1 <- draw.key(list(lines = list( type = c("p",rep("l",length(x))), pch=16, col = c(1,col.vec)), cex = .9, text = list(lab = c("Observations", leg.vec)), draw=TRUE)) }
		if(col.vec[1] == col.vec[2]){
	for( i in 1:length(x)) { 
		panel.points(age.vec, x[[i]][age.vec+1-min(age.vec),yy], type="l", col=col.vec[i], lty = lty.vec[i])
		panel.points(age.vec, dup[[i]][age.vec+1-min(age.vec),yy], type="l", col=gray(.5), lty = lty.vec[i])
		panel.points(age.vec, dlow[[i]][age.vec+1-min(age.vec),yy], type="l", col=gray(.5), lty = lty.vec[i])
		}
	trellis.unfocus()
	k1 <- draw.key(list(lines = list( type = c("p",rep("l",length(x))), pch=16, col = c(1,col.vec), lty = c(NA, lty.vec)), cex = .9, text = list(lab = c("Observations", leg.vec)), draw=TRUE)) }
	pushViewport(viewport(x = unit(k[1], "mm"), y = unit(k[2], "mm"), width = unit(1, "grobwidth", k1), height = unit(1, "grobheight", k1)))
	grid.draw(k1) 
	}

## ------------------------------------------------------------------------ ##
##  Comparison of the evolution of the periodic life expectancy             ##
## ------------------------------------------------------------------------ ##

.ComparisonTrendsMethods = function(ListValidationLevel3, AgeCoh, Pop, ColorComp, LtyComp){
	NamesOfMethods <- ColVec <- LtyVec <- vector(,length(ListValidationLevel3))
	AgeList <- vector("list",length(ListValidationLevel3))
	for(i in 1:length(ListValidationLevel3)){
		NamesOfMethods[i] <- ListValidationLevel3[[i]][[1]][[1]]$NameMethod
		if(NamesOfMethods[i] == "Method1"){
			ColVec[i] <- ColorComp[1]
			LtyVec[i] <- LtyComp[1]
			}
		if(NamesOfMethods[i] == "Method2"){
			ColVec[i] <- ColorComp[2]
			LtyVec[i] <- LtyComp[2]
			}
		if(NamesOfMethods[i] == "Method3"){
			ColVec[i] <- ColorComp[3]
			LtyVec[i] <- LtyComp[3]
			}
		if(NamesOfMethods[i] == "Method4"){
			ColVec[i] <- ColorComp[4]
			LtyVec[i] <- LtyComp[4]
			}
	}
	PLEObserved <- PLEComb <- vector("list",length(ListValidationLevel3))
	names(PLEObserved) <- names(PLEComb) <- NamesOfMethods
	LengthAgeComb <- vector(,length(ListValidationLevel3))
	for(j in 1:length(ListValidationLevel3)){
		LengthAgeComb[j] <- length(rownames(ListValidationLevel3[[j]][["PerLifeExpFitted"]][[Pop]]$PerLifeExp))
		}
	Id <- which(LengthAgeComb == max(LengthAgeComb))[1]
	AgeComp <- as.numeric(rownames(ListValidationLevel3[[Id]][["PerLifeExpFitted"]][[Pop]]$PerLifeExp))
	if(ListValidationLevel3[[1]][["PerLifeExpFitted"]][[Pop]]$NameMethod != "Method4"){		
		PLEObserved <- ListValidationLevel3[[1]][["PerLifeExpObserved"]][[Pop]][AgeCoh-min(AgeComp)+1,]	
		}
		if(ListValidationLevel3[[1]][["PerLifeExpFitted"]][[Pop]]$NameMethod == "Method4"){		
		PLEObserved <- ListValidationLevel3[[1]][["PerLifeExpObserved"]][[Pop]][AgeCoh +1-min(as.numeric(rownames(ListValidationLevel3[[1]][["PerLifeExpObserved"]][[Pop]]))),]	
		}
	for(j in 1:length(ListValidationLevel3)){
		if(ListValidationLevel3[[j]][["PerLifeExpFitted"]][[Pop]]$NameMethod != "Method4"){	
			PLEComb[[j]] <- ListValidationLevel3[[j]][["PerLifeExpComb"]][[Pop]][AgeCoh-min(AgeComp)+1,]
			}
		if(ListValidationLevel3[[j]][["PerLifeExpFitted"]][[Pop]]$NameMethod == "Method4"){
			PLEComb[[j]] <- ListValidationLevel3[[j]][["PerLifeExpComb"]][[Pop]][AgeCoh+1-min(as.numeric(rownames( ListValidationLevel3[[j]][["PerLifeExpComb"]][[Pop]]))),]
			}
		}
	Lim <- c(min(unlist(PLEComb))-.5,max(unlist(PLEComb))+.5)
	print(.PlotPLEMethods(PLEObserved, PLEComb, AgeCoh , c(55,145), ColVec, LtyVec, NamesOfMethods, as.numeric(names(PLEComb[[1]])), Lim, Pop)) 
	}

.PlotPLEMethods = function(x1, x2, aa, k, col.vec, lty.vec, leg.vec, year, Lim, Pop) {
	expr <- paste(Pop,", comparison of the mortality trends for age ", aa, sep="")
	xx <- x2[[1]]
	nn <- length(c(xx)); gr <- 1; grp1 <- as.factor(rep(gr, nn))
	dat <- as.data.frame(cbind(c(xx), grp1))
	zzz <- with(dat, xyplot(c(xx) ~ rep(1:length(x2[[1]]), length(gr)) | as.factor(grp1), 
	data = dat, layout = c(1,1), as.table = T, xlab = 'Year',
	ylab = "Periodic life expectancies", type= "n", scales = list(x = list(alternating = 1, tck = c(1,0),
	at = seq(1,length(x2[[1]]),10), labels = seq(min(year),max(year),10)), y = list(alternating = 3, 
	tck = c(1,1), limits= Lim)), par.settings = list(fontsize = list(text = 15, points = 7)),
	between=list(y = 0.6, x = 0.6), par.strip.text = list(cex = 1, lines = 1.1),
	strip = strip.custom(which.given = 1, factor.levels = expr, bg ="lightgrey")))
	print(zzz)
	trellis.focus("panel", 1, 1, highlight = FALSE)
	if(col.vec[1] != col.vec[2]){
	panel.points(1:length(x1), x1, type="b", col=1, lty = 2, pch=16)
	for(i in 1:length(x2)){
			panel.points(length(x1):length(x2[[1]]), x2[[i]][length(x1):length(x2[[1]])], type="l", col=col.vec[i])
	}
	trellis.unfocus()
	k1 <- draw.key(list(lines = list( type = c("b",rep("l",length(x2))), pch=16, col = c(1,col.vec), lty = c(2,rep(1,length(x2)))), cex = .9, text = list(lab = c("Observations", leg.vec)), draw=TRUE))
	}
	if(col.vec[1] == col.vec[2]){
	panel.points(1:length(x1), x1, type="b", col=gray(.5), lty = 2, pch=16)
	for(i in 1:length(x2)){
			panel.points(length(x1):length(x2[[1]]), x2[[i]][length(x1):length(x2[[1]])], type="l", col=col.vec[i], lty=lty.vec[i])
	}
	trellis.unfocus()
	k1 <- draw.key(list(lines = list( type = c("b",rep("l",length(x2))), pch=16, col = c(gray(.5),col.vec), lty = c(2,lty.vec)), cex = .9, text = list(lab = c("Observations", leg.vec)), draw=TRUE))
	}
	pushViewport(viewport(x = unit(k[1], "mm"), y = unit(k[2], "mm"), width = unit(1, "grobwidth", k1), height = unit(1, "grobheight", k1)))
	grid.draw(k1)
	}
