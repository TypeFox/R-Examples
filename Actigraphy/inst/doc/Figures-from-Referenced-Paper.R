### R code from vignette source 'Figures-from-Referenced-Paper.Rnw'

###################################################
### code chunk number 1: Figures-from-Referenced-Paper.Rnw:14-16
###################################################
	library(Actigraphy)
	library(lattice)


###################################################
### code chunk number 2: Figures-from-Referenced-Paper.Rnw:20-46
###################################################
	### Load Data 
	data(weekday) 
	
	### Data Management 
	data2 <- NULL 
	data2$act <- weekday[,3] 
	data2$t <- weekday[,2] 
	data2$day <- weekday[,1] 
	data2$date <- factor(data2$day, levels=c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday")) 
	
	### Plot Options and Parameters 
	lb <- c("Midnight", "6AM", "Noon", "6PM", "Midnight") 
	L <- 1440 
	xat <- c(0, L/4, L/2, 3*L/4, L) 
	
	### Plot Figure 
	xyplot(act ~ t | date, data=data2, as.table=TRUE, 
		main="Subject 002 Activity from Monday to Friday", 
		scales=list(x=list(at=xat, labels=lb)), cex.main=0.5, 
		layout=c(1, 5, 1), xlim=c(0, L), xlab="(a)", 
		ylab="Activity", panel=function(x,y) {
			fbase <- create.fourier.basis(rangeval=c(0,L), nbasis=9) 
			fpar <- fdPar(fbase) 
			fd <- smooth.basis(c(1:L), y, fpar) 
			panel.xyplot(x, y, type="h") 
		}) 


###################################################
### code chunk number 3: Figures-from-Referenced-Paper.Rnw:51-84
###################################################
	### Load Data 
	data(act_8pt) 
	data(clinic_8pt) 
	
	### Plot Options and Parameters 
	lb <- c("Midnight", "6AM", "Noon", "6PM", "Midnight") 
	L <- 1440 
	xat <- c(0, L/4, L/2, 3*L/4, L) 
	
	matchid <- fda.matchid(act_8pt[,-1], clinic_8pt, type="factor", grouplab=c("AHI", "NO AHI")) 
	idhigh <- paste("Subj", colnames(matchid$mat)[matchid$cov$"NO AHI" != 1]) 
	idlow <- paste("Subj", colnames(matchid$mat)[matchid$cov$"NO AHI" == 1]) 
	idorder <- c(idhigh, idlow) 
	
	datavec <- matchid$mat 
	dim(datavec) <- NULL
	
	datanew <- data.frame(y=datavec, id=rep(paste("Subj", colnames(matchid$mat)), each=L), t=rep(c(1:L), 8)) 
	datanew$id <- factor(datanew$id, idorder) 
	
	### Plot Figure 
	xyplot(y~t|id, data=datanew, as.table=TRUE, 
		main="Circadian Activity from 8 Subjects", 
		ylab="Activity", xlab="", cex.main=.7, 
		scales=list(x=list(at=xat, labels=lb)), cex=.05, 
		type="p", layout=c(4, 2, 1), ylim=c(0, 1200), 
		xlim=c(0, L), panel=function(x,y) { 
			fbase <- create.fourier.basis(rangeval=c(0, L), nbasis=9) 
			fpar <- fdPar(fbase) 
			sm <- smooth.basis(c(1:L), y, fpar) 
			panel.xyplot(x, y, col=1, cex=0.1) 
			panel.lines(predict(sm$fd,c(1:L)), col=2, lwd=3) 
		}) 


###################################################
### code chunk number 4: Figures-from-Referenced-Paper.Rnw:89-155
###################################################
	### Load Data 
	data(act_8pt) 
	data(clinic_8pt) 
	
	ahidatav2 <- fda.matchid(act_8pt[,-1], clinic_8pt, type="factor", grouplab=c("AHI", "NO AHI")) 
	tempv2 <- ahidatav2[[2]] 
	tempv2[,3] <- ifelse(tempv2[,3] == 0, -1, 1) 
	ahidatav2$cov <- data.frame(id=tempv2$id, mean=1, ahi=tempv2[,3]) 
	
	colv2 <- ifelse(tempv2[,3] == -1, 4, 2) 
	smoothDatav2 <- fda.smoothdata(ahidatav2) 
	geftahiv2 <- flm_cate(smoothDatav2) 
	meanefv2 <- geftahiv2$freg$betaestlist[[1]] 
	ahiefv2 <- geftahiv2$freg$betaestlist[[2]] 
	
	### Plot Options and Parameters 
	L <- 1440 
	xat <- c(0, L/4, L/2, 3*L/4, L) 
	lb <- c("Midnight", "6AM", "Noon", "6PM", "Midnight") 
	
	### Plot Figure 
	par(mfrow=c(2,1), mar=c(4,4,3,1)) 
	plot(0, 0, xlim=c(0,L), ylim=c(0,1200), xaxt="n", xlab="(a)", ylab="Acitivity", type="n", main="Circadian Activity Curves of 8 Subjects") 
	
	for(i in 1:8)
		lines(predict(smoothDatav2$fd$fd, c(1:L))[,i], col=colv2[i]) 
	
	### Plot the group mean activities 
	lines(meanefv2$fd-ahiefv2$fd, col=4, lwd=3) 
	lines(meanefv2$fd+ahiefv2$fd, col=2, lwd=3) 
	
	### Plot the overall mean 
	lines(meanefv2$fd, col=1, lwd=3) 
	
	### Add the axis and legend to finish the plot 
	axis(1, at=xat, labels=lb) 
	legend("topleft", c("AHI High Curves", "AHI High Mean", "AHI Low Curves", "AHI Low Mean ", "Overall Mean"), 
			lty=1, col=c(4,4,2,2,1), lwd=c(1,3,1,3,3), cex=.8) 
	
	### F Test 
	cov2 <- smoothDatav2$cov[, -1] 
	grp2 <- ncol(cov2)
	fd <- smoothDatav2$fd 
	L <- length(fd$argvals) 
	npt <- ncol(fd$y)
	
	fbase <- create.fourier.basis(rangeval=c(0, 1440), nbasis=9) 
	fpar <- fdPar(fbase) 
	xfdlist <- vector("list", grp2) 
	xfdlist[[1]] <- cov2[, 1] + 0 
	
	for(i in 2:grp2)
		xfdlist[[i]] <- cov2[, i] + 0 

	betalist <- xfdlist 
	for(i in 1:grp2)
		betalist[[i]] <- fpar 

	freg2 <- fRegress(fd$fd, xfdlist, betalist) 
	preact2 <- predict(freg2$yhatfdobj, c(1:L)) 
	resid2 <- fd$y - preact2[, 1:npt] 
	sigma2 <- cov(t(resid2)) 
	fregstd2 <- fRegress.stderr(freg2, fd$y2cMap, sigma2) 
	
	Fratio <- Ftest(fd$fd, xfdlist, betalist, argvals = c(1:1440),  nperm=10, xaxt="n") 
	axis(1, at=xat, labels=lb) 


