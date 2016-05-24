### R code from vignette source 'vignette.Rnw'

###################################################
### code chunk number 1: egi129EOGcorrection (eval = FALSE)
###################################################
## library(icaOcularCorrection)
## load("data/mc12.eeg.fil.mat.reshp.rda")
## chan <- colnames(egi129)[1:129]
## #
## res <- icac(x = egi129, channel = chan, noise.sig = c("E14",
## 	"E21", "E126", "E127"))
## save(res, file = "models/res.rda", compress = "xz")


###################################################
### code chunk number 2: loadEGI129res (eval = FALSE)
###################################################
## library(icaOcularCorrection)
## load("models/res.rda")
## #
## load("data/mc12.eeg.fil.mat.reshp.rda")
## chan<-colnames(egi129)[1:129]


###################################################
### code chunk number 3: egi129SummaryTopTen (eval = FALSE)
###################################################
## smry  <- summary(res, print = FALSE)
## save(smry, file = "smryEGI129.rda")
## #
## smry  <- summary(res, ic = 52, print = FALSE)
## save(smry, file = "smryEGI129IC52.rda")
## #
## smry  <- summary(res, ic = 6, print = FALSE)
## save(smry, file = "smryEGI129IC6.rda")
## #
## smry  <- summary(res, ic = 43, print = FALSE)
## save(smry, file = "smryEGI129IC43.rda")


###################################################
### code chunk number 4: ShowEGI129SummaryTopTen
###################################################
load("smryEGI129.rda")
smry[1:10,]


###################################################
### code chunk number 5: loadEGI129resIC101
###################################################
load("smryEGI129IC52.rda")
smry


###################################################
### code chunk number 6: plotIC52E14 (eval = FALSE)
###################################################
## # IC 52 mostly correlates with channel E14
## pdf(file="IC52E14.pdf")
## plot_tric(res, ic = 52, noise.sig = "E14",
## 	trials = 85:104, n.win = 21, new.page = FALSE)
## dev.off()


###################################################
### code chunk number 7: plotNIC52E14 (eval = FALSE)
###################################################
## pdf("plotNICIC52E14.pdf")
## par(mfrow = c(3, 3))
## for(i in 85:93){
## 	plot_nic(x = res, data = egi129, ic = 52, 
## 		trial = i, noise.sig = "E14")
## }
## par(mfrow = c(1, 1))
## dev.off()


###################################################
### code chunk number 8: topomapIC52 (eval = FALSE)
###################################################
## pdf("topomapIC52.pdf")
## topo_ic(x = res, ic = 52, coords = "egi.129") 
## dev.off()


###################################################
### code chunk number 9: loadSummaryIC6
###################################################
load("smryEGI129IC6.rda")
smry


###################################################
### code chunk number 10: plotIC6E126 (eval = FALSE)
###################################################
## pdf(file="IC6E126.pdf")
## plot_tric(res, ic = 6, noise.sig = "E126",
## 	trials = 253:273, n.win = 21, new.page = FALSE)
## dev.off()


###################################################
### code chunk number 11: plotNIC6E126 (eval = FALSE)
###################################################
## pdf("plotNICIC6E126.pdf")
## par(mfrow = c(3, 3))
## for(i in 253:261){
## 	plot_nic(x = res, data = egi129, ic = 6, 
## 		trial = i, noise.sig = "E126")
## }
## par(mfrow = c(1, 1))
## dev.off()


###################################################
### code chunk number 12: topomapIC6 (eval = FALSE)
###################################################
## pdf("topomapIC6")
## topo_ic(x = res, ic = 6, coords = "egi.129") 
## dev.off()


###################################################
### code chunk number 13: loadSummaryIC43
###################################################
load("smryEGI129IC43.rda")
smry


###################################################
### code chunk number 14: plotIC43E21 (eval = FALSE)
###################################################
## pdf(file="IC43E21.pdf")
## plot_tric(res, ic = 43, noise.sig = "E21",
## 	trials = 1:21, n.win = 21, new.page = FALSE)
## dev.off()


###################################################
### code chunk number 15: plotNIC43E21 (eval = FALSE)
###################################################
## pdf("plotNICIC43E21.pdf")
## par(mfrow = c(3, 3))
## for(i in 1:9){
## 	plot_nic(x = res, data = egi129, ic = 43, 
## 		trial = i, noise.sig = "E21")
## }
## par(mfrow = c(1, 1))
## dev.off()


###################################################
### code chunk number 16: topomapIC43 (eval = FALSE)
###################################################
## pdf("topomapIC43")
## topo_ic(x = res, ic = 43, coords = "egi.129") 
## dev.off()


###################################################
### code chunk number 17: blinksegi129 (eval = FALSE)
###################################################
## # you'll need library eRp for this.
## library(eRp)
## load("data/mc12.eeg.fil.mat.reshp.rda")
## # get peaks for egi129
## peaks.egi129 <- get.peaks(egi129,"E21",NULL)
## save(peaks.egi129, file = "data/peaks.egi129.rda", compress= "xz")
## #
## # insert event code 777 at each peak
## egi129$EventCode <- as.character(egi129$EventCode)
## pb<-txtProgressBar(min=1,max=length(peaks.egi129),char="=",
## 	style=3)
## for(i in 1:length(peaks.egi129)){
## 	setTxtProgressBar(pb,i)
## 	tmp <- peaks.egi129[[i]]
## 	if(!is.na(tmp[1])){
## 		for(j in 1:length(tmp)){
## 			egi129[egi129$Trial==i & egi129$Time==tmp[j], 
## 				"EventCode"] <- "777"
## 		}
## 	}
## }
## close(pb)
## #
## # save event codes to later merge with corrected data frame
## evts <- egi129$EventCode
## save(evts, file = "data/evts.peaks.egi129.rda", compress = "xz")
## #
## # grab a 200 ms window around each peak and
## # put into data frame
## x <- as.numeric(rownames(egi129[egi129$EventCode=="777",]))
## x1 <- x - 25
## x2 <- x + 25
## x <- cbind(x1, x2)
## tmp <- egi129[x[1, 1]:x[1, 2], ]
## tmp[1,]$EventCode <- "111111"
## tmp[nrow(tmp),]$EventCode <- "222222"
## pb <- txtProgressBar(min = 1, max = nrow(x), char = "=",
## 	style = 3)
## for(i in 2:nrow(x)){
## 	setTxtProgressBar(pb, i)
## 	tmp1 <- egi129[x[i, 1]:x[i, 2], ]
## 	tmp1[1,]$EventCode <- "111111"
## 	tmp1[nrow(tmp1),]$EventCode <- "222222"
## 	tmp <- rbind(tmp, tmp1)
## }
## close(pb)
## #
## # reset time with t = 0 at event code 777
## rownames(tmp) <- 1:nrow(tmp)
## tmp <- add.time2(x = tmp, markers = list(begin = "111111", 
##      ref = "777", finish = "222222"), sampling.rate = 250)
## eog.uncor <- tmp
## save(eog.uncor, file = "data/eog.uncor.icaOC.rda", compress = "xz")
## rm(tmp); gc(TRUE, TRUE)
## #
## # compute blink average for corrected data
## datc <- res$data
## datc$EventCode <- evts
## #
## # grab a 200 ms window around each peak and
## # put into data frame
## x <- as.numeric(rownames(datc[datc$EventCode=="777",]))
## x1 <- x - 25
## x2 <- x + 25
## x <- cbind(x1, x2)
## tmp <- datc[x[1, 1]:x[1, 2], ]
## tmp[1,]$EventCode <- "111111"
## tmp[nrow(tmp),]$EventCode <- "222222"
## pb <- txtProgressBar(min = 1, max = nrow(x), char = "=",
## 	style = 3)
## for(i in 2:nrow(x)){
## 	setTxtProgressBar(pb, i)
## 	tmp1 <- datc[x[i, 1]:x[i, 2], ]
## 	tmp1[1,]$EventCode <- "111111"
## 	tmp1[nrow(tmp1),]$EventCode <- "222222"
## 	tmp <- rbind(tmp, tmp1)
## }
## close(pb)
## #
## # reset time with t = 0 at event code 777
## rownames(tmp) <- 1:nrow(tmp)
## tmp <- add.time2(x = tmp, markers = list(begin = "111111", 
##      ref = "777", finish = "222222"), sampling.rate = 250)
## eog.cor <- tmp
## save(eog.cor, file = "data/eog.cor.icaOC.rda", compress = "xz")
## rm(tmp); gc(TRUE, TRUE)
## #
## # topomap for uncorrected
## avg <- tapply(eog.uncor[,1], eog.uncor$Time, mean)
## chan<-colnames(egi129)[1:129]
## avg.dat<-data.frame(Time=as.numeric(names(avg)),Amplitude=avg,Channel=chan[1])
## pb <- txtProgressBar(min = 1, max = length(chan)-1, char = "=",
## 	style = 3)
## for(i in 2:length(chan)){
## 	setTxtProgressBar(pb, i)
## 	avg <- tapply(eog.uncor[,chan[i]], eog.uncor$Time, mean)
## 	avg.dat<-rbind(avg.dat,data.frame(Time=as.numeric(names(avg)),
## 		Amplitude=avg,Channel=chan[i]))
## }
## close(pb)
## avg.dat$Channel<-as.factor(avg.dat$Channel)
## coords<-des("egi.129")$cart
## avg.dat.uncor<-merge(avg.dat,coords[,c("x","y","Channel")],by="Channel")
## m.uncor<-gam(Amplitude~te(x,y,bs="ts",k=11),dat=avg.dat.uncor)
## #
## # topomap for corrected
## avg <- tapply(eog.cor[,1], eog.cor$Time, mean)
## chan<-colnames(egi129)[1:129]
## avg.dat<-data.frame(Time=as.numeric(names(avg)),Amplitude=avg,Channel=chan[1])
## pb <- txtProgressBar(min = 1, max = length(chan)-1, char = "=",
## 	style = 3)
## for(i in 2:length(chan)){
## 	setTxtProgressBar(pb, i)
## 	avg <- tapply(eog.cor[,chan[i]], eog.cor$Time, mean)
## 	avg.dat<-rbind(avg.dat,data.frame(Time=as.numeric(names(avg)),
## 		Amplitude=avg,Channel=chan[i]))
## }
## close(pb)
## avg.dat$Channel<-as.factor(avg.dat$Channel)
## coords<-des("egi.129")$cart
## avg.dat.cor<-merge(avg.dat,coords[,c("x","y","Channel")],by="Channel")
## m.cor<-gam(Amplitude~te(x,y,bs="ts",k=11),dat=avg.dat.cor)
## #
## # get plotting info
## pi.uncor<-plotGAM(m.uncor,too.far=des("egi.129")$too.far,plot=FALSE)
## pi.cor<-plotGAM(m.cor,too.far=des("egi.129")$too.far,plot=FALSE)
## #
## # waveforms
## avg.uncor <- tapply(eog.uncor$E21, eog.uncor$Time, mean)
## avg.cor <- tapply(eog.cor$E21, eog.cor$Time, mean)
## time <- as.numeric(names(avg.uncor))
## #
## # create plot
## pdf(file = "uncorrectedCorrected0.pdf")
## par(mfrow = c(2, 2), mar = c(5.1, 4.1, 4.1, 2.1))
## # plot waveforms
## plot(time, avg.uncor - min(avg.uncor), type = "l",
## 	xlab = "Time (ms)", ylab = "Amplitude",ylim=c(0,400),
## 	main = "E21")
## lines(time, avg.cor - min(avg.cor), col = 2)
## legend("topleft", legend = c("blinks in uncorrected data",
## 	"blinks in corrected data"), lty = 1, col = 1:2, bty = "n")
## par(mar = c(1, 1, 1, 1))
## zlimit <- range(rbind(pi.uncor$mat, pi.cor$mat), na.rm = TRUE)
## # skip a plotting region
## plot.new()
## # plot uncorrected topomap
## image(pi.uncor$xm, pi.uncor$xm, pi.uncor$mat, col = topo.colors(100),
## 	xlab = "", ylab = "", zlim = zlimit, axes = FALSE)
## contour(pi.uncor$xm, pi.uncor$xm, pi.uncor$mat, add = TRUE)
## title(main = "Uncorrected", line = -12.5)
## # plot corrected topomap
## image(pi.cor$xm, pi.cor$xm, pi.cor$mat, col = topo.colors(100),
## 	xlab = "", ylab = "", zlim = zlimit, axes = FALSE)
## contour(pi.cor$xm, pi.cor$xm, pi.cor$mat, add = TRUE)
## title(main = "Corrected", line = -12.5)
## par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
## dev.off()


###################################################
### code chunk number 18: lookAtNumTrials3 (eval = FALSE)
###################################################
## pdf(file="ICNumTrials3.pdf")
## plot(smry$NumTrial, type = "h", xlab="IC", ylab = "Number of Trials", 
## 	xaxt = "n")
## myat<-pretty(1:nrow(smry), 10)
## myat[1]<-1
## myat<-myat[1:(length(myat)-1)]
## mylab<-smry$IC[myat]
## axis(side = 1, at = myat, labels = mylab, cex.axis = 0.85)
## abline(v = 3, col = 2, lty = 3)
## text(x = 4, y = 1000, labels = paste("IC",smry$IC[3]), cex = 0.85, adj = 0)
## dev.off()


###################################################
### code chunk number 19: zerooutFirstBreakPoint (eval = FALSE)
###################################################
## my.what <- list()
## for(i in 1:3){
## 	my.what[[i]] <- c(smry$IC[i], "-")
## }
## res.up <- update(object = res, what = my.what)
## save(res.up, file = "models/res.up.rda", compress = "xz")


###################################################
### code chunk number 20: lookAtBlinksAgain (eval = FALSE)
###################################################
## # compute blink average for corrected data
## datc <- res.up$data
## datc$EventCode <- evts
## #
## # grab a 200 ms window around each peak and
## # put into data frame
## x <- as.numeric(rownames(datc[datc$EventCode=="777",]))
## x1 <- x - 25
## x2 <- x + 25
## x <- cbind(x1, x2)
## tmp <- datc[x[1, 1]:x[1, 2], ]
## tmp[1,]$EventCode <- "111111"
## tmp[nrow(tmp),]$EventCode <- "222222"
## pb <- txtProgressBar(min = 1, max = nrow(x), char = "=",
## 	style = 3)
## for(i in 2:nrow(x)){
## 	setTxtProgressBar(pb, i)
## 	tmp1 <- datc[x[i, 1]:x[i, 2], ]
## 	tmp1[1,]$EventCode <- "111111"
## 	tmp1[nrow(tmp1),]$EventCode <- "222222"
## 	tmp <- rbind(tmp, tmp1)
## }
## close(pb)
## #
## # reset time with t = 0 at event code 777
## rownames(tmp) <- 1:nrow(tmp)
## tmp <- add.time2(x = tmp, markers = list(begin = "111111", 
##      ref = "777", finish = "222222"), sampling.rate = 250)
## eog.cor2 <- tmp
## save(eog.cor2, file = "data/eog.cor2.icaOC.rda", compress = "xz")
## rm(tmp); gc(TRUE, TRUE)
## #
## # topomap for corrected
## avg <- tapply(eog.cor2[,1], eog.cor2$Time, mean)
## chan<-colnames(egi129)[1:129]
## avg.dat<-data.frame(Time=as.numeric(names(avg)),Amplitude=avg,Channel=chan[1])
## pb <- txtProgressBar(min = 1, max = length(chan)-1, char = "=",
## 	style = 3)
## for(i in 2:length(chan)){
## 	setTxtProgressBar(pb, i)
## 	avg <- tapply(eog.cor2[,chan[i]], eog.cor2$Time, mean)
## 	avg.dat<-rbind(avg.dat,data.frame(Time=as.numeric(names(avg)),
## 		Amplitude=avg,Channel=chan[i]))
## }
## close(pb)
## avg.dat$Channel<-as.factor(avg.dat$Channel)
## coords<-des("egi.129")$cart
## avg.dat.cor2<-merge(avg.dat,coords[,c("x","y","Channel")],by="Channel")
## m.cor2<-gam(Amplitude~te(x,y,bs="ts",k=11),dat=avg.dat.cor2)
## #
## # get plotting info
## pi.cor2<-plotGAM(m.cor2,too.far=des("egi.129")$too.far,plot=FALSE)
## #
## # waveforms
## avg.cor2 <- tapply(eog.cor2$E21, eog.cor2$Time, mean)
## time <- as.numeric(names(avg.uncor))
## #
## # create plot
## pdf(file = "uncorrectedCorrected3.pdf")
## par(mfrow = c(2, 2), mar = c(5.1, 4.1, 4.1, 2.1))
## # plot waveforms
## plot(time, avg.uncor - min(avg.uncor), type = "l",
## 	xlab = "Time (ms)", ylab = "Amplitude",ylim=c(0,400),
## 	main = "E21")
## lines(time, avg.cor2 - min(avg.cor2), col = 2)
## legend("topleft", legend = c("blinks in uncorrected data",
## 	"blinks in corrected data"), lty = 1, col = 1:2, bty = "n")
## par(mar = c(1, 1, 1, 1))
## zlimit <- range(rbind(pi.uncor$mat, pi.cor$mat), na.rm = TRUE)
## # skip a plotting region
## plot.new()
## # plot uncorrected topomap
## image(pi.uncor$xm, pi.uncor$xm, pi.uncor$mat, col = topo.colors(100),
## 	xlab = "", ylab = "", zlim = zlimit, axes = FALSE)
## contour(pi.uncor$xm, pi.uncor$xm, pi.uncor$mat, add = TRUE)
## title(main = "Uncorrected", line = -12.5)
## # plot corrected topomap
## image(pi.cor2$xm, pi.cor2$xm, pi.cor2$mat, col = topo.colors(100),
## 	xlab = "", ylab = "", zlim = zlimit, axes = FALSE)
## contour(pi.cor2$xm, pi.cor2$xm, pi.cor2$mat, add = TRUE)
## title(main = "Corrected", line = -12.5)
## par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
## dev.off()


###################################################
### code chunk number 21: lookAtNumTrials5 (eval = FALSE)
###################################################
## pdf(file="ICNumTrials5.pdf")
## plot(smry$NumTrial, type = "h", xlab="IC", ylab = "Number of Trials", 
## 	xaxt = "n")
## myat<-pretty(1:nrow(smry), 10)
## myat[1]<-1
## myat<-myat[1:(length(myat)-1)]
## mylab<-smry$IC[myat]
## axis(side = 1, at = myat, labels = mylab, cex.axis = 0.85)
## abline(v = 5, col = 2, lty = 3)
## text(x = 6, y = 1000, labels = paste("IC",smry$IC[5]), cex = 0.85, adj = 0)
## dev.off()


###################################################
### code chunk number 22: zerooutSecondBreakPoint (eval = FALSE)
###################################################
## my.what <- list()
## for(i in 1:5){
## 	my.what[[i]] <- c(smry$IC[i], "-")
## }
## res.up2 <- update(object = res, what = my.what)
## save(res.up2, file = "models/res.up2.rda", compress = "xz")


###################################################
### code chunk number 23: lookAtBlinksAgain (eval = FALSE)
###################################################
## # compute blink average for corrected data
## datc <- res.up$data
## datc$EventCode <- evts
## #
## # grab a 200 ms window around each peak and
## # put into data frame
## x <- as.numeric(rownames(datc[datc$EventCode=="777",]))
## x1 <- x - 25
## x2 <- x + 25
## x <- cbind(x1, x2)
## tmp <- datc[x[1, 1]:x[1, 2], ]
## tmp[1,]$EventCode <- "111111"
## tmp[nrow(tmp),]$EventCode <- "222222"
## pb <- txtProgressBar(min = 1, max = nrow(x), char = "=",
## 	style = 3)
## for(i in 2:nrow(x)){
## 	setTxtProgressBar(pb, i)
## 	tmp1 <- datc[x[i, 1]:x[i, 2], ]
## 	tmp1[1,]$EventCode <- "111111"
## 	tmp1[nrow(tmp1),]$EventCode <- "222222"
## 	tmp <- rbind(tmp, tmp1)
## }
## close(pb)
## #
## # reset time with t = 0 at event code 777
## rownames(tmp) <- 1:nrow(tmp)
## tmp <- add.time2(x = tmp, markers = list(begin = "111111", 
##      ref = "777", finish = "222222"), sampling.rate = 250)
## eog.cor2 <- tmp
## save(eog.cor2, file = "data/eog.cor2.icaOC.rda", compress = "xz")
## rm(tmp); gc(TRUE, TRUE)
## #
## # topomap for corrected
## avg <- tapply(eog.cor2[,1], eog.cor2$Time, mean)
## chan<-colnames(egi129)[1:129]
## avg.dat<-data.frame(Time=as.numeric(names(avg)),Amplitude=avg,Channel=chan[1])
## pb <- txtProgressBar(min = 1, max = length(chan)-1, char = "=",
## 	style = 3)
## for(i in 2:length(chan)){
## 	setTxtProgressBar(pb, i)
## 	avg <- tapply(eog.cor2[,chan[i]], eog.cor2$Time, mean)
## 	avg.dat<-rbind(avg.dat,data.frame(Time=as.numeric(names(avg)),
## 		Amplitude=avg,Channel=chan[i]))
## }
## close(pb)
## avg.dat$Channel<-as.factor(avg.dat$Channel)
## coords<-des("egi.129")$cart
## avg.dat<-merge(avg.dat,coords[,c("x","y","Channel")],by="Channel")
## m.cor2<-gam(Amplitude~te(x,y,k=11),dat=avg.dat)
## #
## # get plotting info
## pi.cor2<-plotGAM(m.cor2,too.far=des("egi.129")$too.far,plot=FALSE)
## #
## # waveforms
## avg.cor2 <- tapply(eog.cor2$E21, eog.cor2$Time, mean)
## time <- as.numeric(names(avg.uncor))
## #
## # create plot
## pdf(file = "uncorrectedCorrected2.pdf")
## par(mfrow = c(2, 2), mar = c(5.1, 4.1, 4.1, 2.1))
## # plot waveforms
## plot(time, avg.uncor - min(avg.uncor), type = "l",
## 	xlab = "Time (ms)", ylab = "Amplitude",ylim=c(0,400),
## 	main = "E21")
## lines(time, avg.cor - min(avg.cor), col = 2)
## legend("topleft", legend = c("blinks in uncorrected data",
## 	"blinks in corrected data"), lty = 1, col = 1:2, bty = "n")
## par(mar = c(1, 1, 1, 1))
## zlimit <- range(rbind(pi.uncor$mat, pi.cor$mat), na.rm = TRUE)
## # skip a plotting region
## plot.new()
## # plot uncorrected topomap
## image(pi.uncor$xm, pi.uncor$xm, pi.uncor$mat, col = topo.colors(100),
## 	xlab = "", ylab = "", zlim = zlimit, axes = FALSE)
## contour(pi.uncor$xm, pi.uncor$xm, pi.uncor$mat, add = TRUE)
## title(main = "Uncorrected", line = -12.5)
## # plot corrected topomap
## image(pi.cor2$xm, pi.cor2$xm, pi.cor2$mat, col = topo.colors(100),
## 	xlab = "", ylab = "", zlim = zlimit, axes = FALSE)
## contour(pi.cor2$xm, pi.cor2$xm, pi.cor2$mat, add = TRUE)
## title(main = "Corrected", line = -12.5)
## par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
## dev.off()


###################################################
### code chunk number 24: lookAtBlinksAgain2 (eval = FALSE)
###################################################
## # compute blink average for corrected data
## datc <- res.up2$data
## datc$EventCode <- evts
## #
## # grab a 200 ms window around each peak and
## # put into data frame
## x <- as.numeric(rownames(datc[datc$EventCode=="777",]))
## x1 <- x - 25
## x2 <- x + 25
## x <- cbind(x1, x2)
## tmp <- datc[x[1, 1]:x[1, 2], ]
## tmp[1,]$EventCode <- "111111"
## tmp[nrow(tmp),]$EventCode <- "222222"
## pb <- txtProgressBar(min = 1, max = nrow(x), char = "=",
## 	style = 3)
## for(i in 2:nrow(x)){
## 	setTxtProgressBar(pb, i)
## 	tmp1 <- datc[x[i, 1]:x[i, 2], ]
## 	tmp1[1,]$EventCode <- "111111"
## 	tmp1[nrow(tmp1),]$EventCode <- "222222"
## 	tmp <- rbind(tmp, tmp1)
## }
## close(pb)
## #
## # reset time with t = 0 at event code 777
## rownames(tmp) <- 1:nrow(tmp)
## tmp <- add.time2(x = tmp, markers = list(begin = "111111", 
##      ref = "777", finish = "222222"), sampling.rate = 250)
## eog.cor3 <- tmp
## save(eog.cor3, file = "data/eog.cor3.icaOC.rda", compress = "xz")
## rm(tmp); gc(TRUE, TRUE)
## #
## # topomap for corrected
## avg <- tapply(eog.cor3[,1], eog.cor3$Time, mean)
## chan<-colnames(egi129)[1:129]
## avg.dat<-data.frame(Time=as.numeric(names(avg)),Amplitude=avg,Channel=chan[1])
## pb <- txtProgressBar(min = 1, max = length(chan)-1, char = "=",
## 	style = 3)
## for(i in 2:length(chan)){
## 	setTxtProgressBar(pb, i)
## 	avg <- tapply(eog.cor3[,chan[i]], eog.cor3$Time, mean)
## 	avg.dat<-rbind(avg.dat,data.frame(Time=as.numeric(names(avg)),
## 		Amplitude=avg,Channel=chan[i]))
## }
## close(pb)
## avg.dat$Channel<-as.factor(avg.dat$Channel)
## coords<-des("egi.129")$cart
## avg.dat.cor3<-merge(avg.dat,coords[,c("x","y","Channel")],by="Channel")
## m.cor3<-gam(Amplitude~te(x,y,bs="ts",k=11),dat=avg.dat.cor3)
## #
## # get plotting info
## pi.cor3<-plotGAM(m.cor3,too.far=des("egi.129")$too.far,plot=FALSE)
## #
## # waveforms
## avg.cor3 <- tapply(eog.cor3$E21, eog.cor3$Time, mean)
## time <- as.numeric(names(avg.uncor))
## #
## # create plot
## pdf(file = "uncorrectedCorrected5.pdf")
## par(mfrow = c(2, 2), mar = c(5.1, 4.1, 4.1, 2.1))
## # plot waveforms
## plot(time, avg.uncor - min(avg.uncor), type = "l",
## 	xlab = "Time (ms)", ylab = "Amplitude",ylim=c(0,400),
## 	main = "E21")
## lines(time, avg.cor3 - min(avg.cor3), col = 2)
## legend("topleft", legend = c("blinks in uncorrected data",
## 	"blinks in corrected data"), lty = 1, col = 1:2, bty = "n")
## par(mar = c(1, 1, 1, 1))
## zlimit <- range(rbind(pi.uncor$mat, pi.cor$mat), na.rm = TRUE)
## # skip a plotting region
## plot.new()
## # plot uncorrected topomap
## image(pi.uncor$xm, pi.uncor$xm, pi.uncor$mat, col = topo.colors(100),
## 	xlab = "", ylab = "", zlim = zlimit, axes = FALSE)
## contour(pi.uncor$xm, pi.uncor$xm, pi.uncor$mat, add = TRUE)
## title(main = "Uncorrected", line = -12.5)
## # plot corrected topomap
## image(pi.cor3$xm, pi.cor3$xm, pi.cor3$mat, col = topo.colors(100),
## 	xlab = "", ylab = "", zlim = zlimit, axes = FALSE)
## contour(pi.cor3$xm, pi.cor3$xm, pi.cor3$mat, add = TRUE)
## title(main = "Corrected", line = -12.5)
## par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
## dev.off()


###################################################
### code chunk number 25: lookAtNumTrials12 (eval = FALSE)
###################################################
## pdf(file="ICNumTrials12.pdf")
## plot(smry$NumTrial, type = "h", xlab="IC", ylab = "Number of Trials", 
## 	xaxt = "n")
## myat<-pretty(1:nrow(smry), 10)
## myat[1]<-1
## myat<-myat[1:(length(myat)-1)]
## mylab<-smry$IC[myat]
## axis(side = 1, at = myat, labels = mylab, cex.axis = 0.85)
## abline(v = 12, col = 2, lty = 3)
## text(x = 13, y = 1000, labels = paste("IC",smry$IC[12]), cex = 0.85, adj = 0)
## dev.off()


###################################################
### code chunk number 26: zerooutThirdBreakPoint (eval = FALSE)
###################################################
## my.what <- list()
## for(i in 1:12){
## 	my.what[[i]] <- c(smry$IC[i], "-")
## }
## res.up3 <- update(object = res, what = my.what)
## save(res.up3, file = "models/res.up3.rda", compress = "xz")


###################################################
### code chunk number 27: lookAtBlinksAgain3 (eval = FALSE)
###################################################
## # compute blink average for corrected data
## datc <- res.up3$data
## datc$EventCode <- evts
## #
## # grab a 200 ms window around each peak and
## # put into data frame
## x <- as.numeric(rownames(datc[datc$EventCode=="777",]))
## x1 <- x - 25
## x2 <- x + 25
## x <- cbind(x1, x2)
## tmp <- datc[x[1, 1]:x[1, 2], ]
## tmp[1,]$EventCode <- "111111"
## tmp[nrow(tmp),]$EventCode <- "222222"
## pb <- txtProgressBar(min = 1, max = nrow(x), char = "=",
## 	style = 3)
## for(i in 2:nrow(x)){
## 	setTxtProgressBar(pb, i)
## 	tmp1 <- datc[x[i, 1]:x[i, 2], ]
## 	tmp1[1,]$EventCode <- "111111"
## 	tmp1[nrow(tmp1),]$EventCode <- "222222"
## 	tmp <- rbind(tmp, tmp1)
## }
## close(pb)
## #
## # reset time with t = 0 at event code 777
## rownames(tmp) <- 1:nrow(tmp)
## tmp <- add.time2(x = tmp, markers = list(begin = "111111", 
##      ref = "777", finish = "222222"), sampling.rate = 250)
## eog.cor4 <- tmp
## save(eog.cor4, file = "data/eog.cor4.icaOC.rda", compress = "xz")
## rm(tmp); gc(TRUE, TRUE)
## #
## # topomap for corrected
## avg <- tapply(eog.cor4[,1], eog.cor4$Time, mean)
## chan<-colnames(egi129)[1:129]
## avg.dat<-data.frame(Time=as.numeric(names(avg)),Amplitude=avg,Channel=chan[1])
## pb <- txtProgressBar(min = 1, max = length(chan)-1, char = "=",
## 	style = 3)
## for(i in 2:length(chan)){
## 	setTxtProgressBar(pb, i)
## 	avg <- tapply(eog.cor4[,chan[i]], eog.cor4$Time, mean)
## 	avg.dat<-rbind(avg.dat,data.frame(Time=as.numeric(names(avg)),
## 		Amplitude=avg,Channel=chan[i]))
## }
## close(pb)
## avg.dat$Channel<-as.factor(avg.dat$Channel)
## coords<-des("egi.129")$cart
## avg.dat.cor4<-merge(avg.dat,coords[,c("x","y","Channel")],by="Channel")
## m.cor4<-gam(Amplitude~te(x,y,bs="ts",k=11),dat=avg.dat.cor4)
## #
## # get plotting info
## pi.cor4<-plotGAM(m.cor4,too.far=des("egi.129")$too.far,plot=FALSE)
## #
## # waveforms
## avg.cor4 <- tapply(eog.cor4$E21, eog.cor4$Time, mean)
## time <- as.numeric(names(avg.uncor))
## #
## # create plot
## pdf(file = "uncorrectedCorrected12.pdf")
## par(mfrow = c(2, 2), mar = c(5.1, 4.1, 4.1, 2.1))
## # plot waveforms
## plot(time, avg.uncor - min(avg.uncor), type = "l",
## 	xlab = "Time (ms)", ylab = "Amplitude",ylim=c(0,400),
## 	main = "E21")
## lines(time, avg.cor4 - min(avg.cor4), col = 2)
## legend("topleft", legend = c("blinks in uncorrected data",
## 	"blinks in corrected data"), lty = 1, col = 1:2, bty = "n")
## par(mar = c(1, 1, 1, 1))
## zlimit <- range(rbind(pi.uncor$mat, pi.cor$mat), na.rm = TRUE)
## # skip a plotting region
## plot.new()
## # plot uncorrected topomap
## image(pi.uncor$xm, pi.uncor$xm, pi.uncor$mat, col = topo.colors(100),
## 	xlab = "", ylab = "", zlim = zlimit, axes = FALSE)
## contour(pi.uncor$xm, pi.uncor$xm, pi.uncor$mat, add = TRUE)
## title(main = "Uncorrected", line = -12.5)
## # plot corrected topomap
## image(pi.cor4$xm, pi.cor4$xm, pi.cor4$mat, col = topo.colors(100),
## 	xlab = "", ylab = "", zlim = zlimit, axes = FALSE)
## contour(pi.cor4$xm, pi.cor4$xm, pi.cor4$mat, add = TRUE)
## title(main = "Corrected", line = -12.5)
## par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
## dev.off()
## #
## # save info
## save(avg.cor,avg.cor2,avg.cor3,avg.cor4,avg.dat.cor,avg.dat.cor2,
## 	avg.dat.cor3,avg.dat.cor4,avg.dat.uncor,avg.uncor,m.cor,m.cor2,
## 	m.cor3,m.cor4,m.uncor,pi.cor,pi.cor2,pi.cor3,pi.cor4,pi.uncor,
## 	file="models/updating_and_topomaps.rda",compress="xz")


###################################################
### code chunk number 28: plotBeforeAfterAverages (eval = FALSE)
###################################################
## # create egi129 net mask
## mat <- pi.cor$mat
## mat[!is.na(mat)]<-1
## rownames(mat) <- pi.cor$xm
## colnames(mat) <- pi.cor$ym
## save(mat, file = "data/egi129NetMask.rda", compress = "xz")
## #
## # plot before and after averages
## pdf(file = "beforeAfter5IC.pdf")
## plot_avgba(x = res.up2, data = egi129, channel = c("E17", "E11",
## 	"E16", "E6", "E129", "E55", "E62", "E72", "E75", "E81"),
## 	ylim = c(-5, 5), new.page = FALSE)
## #
## # add net and channel names
## par(new = TRUE)
## split.screen(c(3,3))
## screen(16)
## par(mar=c(0.5,0.5,0.5,0.5))
## image(as.numeric(rownames(mat)), as.numeric(colnames(mat)), 
## 	mat,col=rgb(190,190,190,100,maxColorValue=255), ann = FALSE,
## 	axes = FALSE)
## coords<-des("egi.129")$cart
## for(i in 1:nrow(coords)){
## 	if(coords[i,"Channel"]%in%c("E17", "E11",
## 	"E16", "E6", "E129", "E55", "E62", "E72", 
## 	"E75", "E81")){
## 		mycol <- "red"
## 	}else{
## 		mycol <- "black"
## 	}
## 	text(coords[i,"x"], coords[i,"y"], coords[i, "Channel"], cex = 0.35,
## 		col = mycol)
## }
## dev.off()


###################################################
### code chunk number 29: formatEEGLABbeforeData (eval = FALSE)
###################################################
## # load data
## dat <- read.table("data/mc12artrejDATA.txt", 
## 	header = TRUE, sep = "\t", stringsAsFactor = FALSE)				  
## chans <- dat[,1]
## chans[length(chans)] <- "E129"
## # transpose so that channels are now columns
## dat <- t(dat)
## dat <- dat[2:nrow(dat),]
## rownames(dat) <- 1:nrow(dat)
## #
## # change amplitudes from character to numeric
## tmp <- cbind(as.numeric(as.character(dat[,1])), 
## 	as.numeric(as.character(dat[,2])))
## for(i in 3:ncol(dat)){
## 	tmp <- cbind(tmp, as.numeric(as.character(dat[,i])))
## }
## #
## # add column names = channel names
## colnames(tmp) <- chans
## rownames(tmp) <- 1:nrow(tmp)
## #
## # coerce to data frame
## dat <- as.data.frame(tmp, stringsAsFactors = FALSE)
## #
## # add time column
## time <- read.table("data/mc12artrejDATA.txt", 
## 	header = FALSE, sep = "\t", stringsAsFactor = FALSE, nrow = 1)
## time <- as.numeric(time[1, ])
## dat$Time <- time[2:length(time)]
## #
## # get rid of NAs
## dat <- na.omit(dat)
## # 
## # add event codes to create epochs
## dat$EventCode <- 0
## dat[dat$Time == -200, "EventCode"] <- 1
## dat[dat$Time == 0, "EventCode"] <- 2
## dat$X <- c(dat$Time[1:(nrow(dat)-1)] - dat$Time[2:(nrow(dat))], -4)
## dat[dat$X > -4, "EventCode"] <- 3
## dat[nrow(dat), "EventCode"] <- 3
## dat <- dat[,1:43]
## dat$EventCode <- as.character(dat$EventCode)
## #
## # add trial
## markers = list(begin = "1", finish = "3")
## dat$EventCode <- as.numeric(as.character(dat$EventCode))
## mytrial <- NULL
## expand.from.to <- function(iii){
## 	return(invisible(iii[1]:iii[2]))
## }
## begin <- which(dat$EventCode == markers$begin)
## finish <- vector("numeric")
## for(jjj in 1:length(begin)){
## 	cat(paste("fetching finishing row number(", jjj, 
## 		" of ", length(begin), ")", sep = ""), "\n")
## 	tmp <- dat[begin[jjj] : nrow(dat), ]$EventCode
## 	names(tmp) <- rownames(dat[begin[jjj] : nrow(dat), ])
## 	finish <- c(finish, as.numeric(names(which(tmp == markers$finish)[1])))
## }
## if(length(begin) != length(finish)){
## 	print(c(table(dat$eeg[dat$eeg[, "EventCode"] == markers$begin, 
## 		"EventCode"]), table(dat$eeg[dat$eeg[, "EventCode"] == markers$finish, 
## 		"EventCode"])))
## 	stop("problem in trigger codes: length(begin) != length(finish)\n")
## }
## dat$Trial <- NA
## tmp <- cbind(begin, finish)
## segment <- list()
## for(kk in 1:nrow(tmp)){
## 	segment[[kk]]<-expand.from.to(tmp[kk,])
## }
## for(i in 1:length(segment)){
## 	cat("adding trial", i, "of", length(segment), "...\n")
## 	dat$Trial[segment[[i]]] <- i
## }
## save(dat, file = "data/mc12artrejDATA.rda",
## 	compress = "xz")


###################################################
### code chunk number 30: formatEEGLABcorrectedData (eval = FALSE)
###################################################
## # load data
## dat <- read.table("data/mc12postICAallchansDATA.txt", 
## 	header = TRUE, sep = "\t", stringsAsFactor = FALSE)				  
## #
## # transpose so that channels are now columns
## dat <- t(dat)
## dat <- dat[2:nrow(dat),]
## rownames(dat) <- 1:nrow(dat)
## #
## # change amplitudes from character to numeric
## tmp <- cbind(as.numeric(as.character(dat[,1])), 
## 	as.numeric(as.character(dat[,2])))
## for(i in 3:ncol(dat)){
## 	tmp <- cbind(tmp, as.numeric(as.character(dat[,i])))
## }
## #
## # add column names = channel names
## colnames(tmp) <- paste("E", 1:129, sep = "")
## rownames(tmp) <- 1:nrow(tmp)
## #
## # coerce to data frame
## dat <- as.data.frame(tmp, stringsAsFactors = FALSE)
## #
## # add time column
## time <- read.table("data/mc12postICAallchansDATA.txt", 
## 	header = FALSE, sep = "\t", stringsAsFactor = FALSE, nrow = 1)
## time <- as.numeric(time[1, ])
## dat$Time <- time[2:length(time)]
## #
## # get rid of NAs
## dat <- na.omit(dat)
## # 
## # add event codes to create epochs
## dat$EventCode <- 0
## dat[dat$Time == -200, "EventCode"] <- 1
## dat[dat$Time == 0, "EventCode"] <- 2
## dat$X <- c(dat$Time[1:(nrow(dat)-1)] - dat$Time[2:(nrow(dat))], -4)
## dat[dat$X > -4, "EventCode"] <- 3
## dat[nrow(dat), "EventCode"] <- 3
## dat <- dat[,1:131]
## dat$EventCode <- as.character(dat$EventCode)
## #
## # add trial
## markers = list(begin = "1", finish = "3")
## dat$EventCode <- as.numeric(as.character(dat$EventCode))
## mytrial <- NULL
## expand.from.to <- function(iii){
## 	return(invisible(iii[1]:iii[2]))
## }
## begin <- which(dat$EventCode == markers$begin)
## finish <- vector("numeric")
## for(jjj in 1:length(begin)){
## 	cat(paste("fetching finishing row number(", jjj, 
## 		" of ", length(begin), ")", sep = ""), "\n")
## 	tmp <- dat[begin[jjj] : nrow(dat), ]$EventCode
## 	names(tmp) <- rownames(dat[begin[jjj] : nrow(dat), ])
## 	finish <- c(finish, as.numeric(names(which(tmp == markers$finish)[1])))
## }
## if(length(begin) != length(finish)){
## 	print(c(table(dat$eeg[dat$eeg[, "EventCode"] == markers$begin, 
## 		"EventCode"]), table(dat$eeg[dat$eeg[, "EventCode"] == markers$finish, 
## 		"EventCode"])))
## 	stop("problem in trigger codes: length(begin) != length(finish)\n")
## }
## dat$Trial <- NA
## tmp <- cbind(begin, finish)
## segment <- list()
## for(kk in 1:nrow(tmp)){
## 	segment[[kk]]<-expand.from.to(tmp[kk,])
## }
## for(i in 1:length(segment)){
## 	cat("adding trial", i, "of", length(segment), "...\n")
## 	dat$Trial[segment[[i]]] <- i
## }
## save(dat, file = "data/mc12postICAallchansDATA.rda",
## 	compress = "xz")


###################################################
### code chunk number 31: EEGLABbeforeAfter (eval = FALSE)
###################################################
## library(eRp)
## library(icaOcularCorrection)
## #
## ###############################################
## # look at data corrected with EEGLAB
## load("data/mc12postICAallchansEEGLAB.rda")
## dat <- dat[dat$Time >= -200 & dat$Time <= 1500, ]
## rownames(dat)<-1:nrow(dat)
## eeglab <- list()
## eeglab$data <- dat
## rm(dat)
## gc(TRUE,TRUE)
## class(eeglab) <- "icac"
## #
## load("data/mc12artrejEEGLAB.rda")
## dat <- dat[dat$Time >= -200 & dat$Time <= 1500, ]
## rownames(dat)<-1:nrow(dat)
## # 
## pdf(file="eeglab.pdf")
## plot_avgba(x = eeglab, data = dat, channel = c("E17", "E11",
## 	"E16", "E6", "E129", "E55", "E62", "E72", "E75", "E81"), 
## 	ylim = c(-5, 5), new.page = FALSE)
## dev.off()
## 
## dev.new()
## plot_trba(x = eeglab, data = dat, channel = "E21")
## 
## cor(eeglab$data$E21, dat$E21)
## # [1] 0.03127444
## #
## pdf(file="eeglab.pdf")
## plot_trba(x = eeglab, data = dat, channel = "E21", new.page = FALSE)
## dev.off()


###################################################
### code chunk number 32: blinksEEGLAB (eval = FALSE)
###################################################
## load("data/mc12postICAallchansEEGLAB.rda")
## dat <- dat[dat$Time >= -200 & dat$Time <= 1500, ]
## rownames(dat) <- 1:nrow(dat)
## eeglab <- list()
## eeglab$data <- dat
## rm(dat)
## gc(TRUE,TRUE)
## class(eeglab) <- "icac"
## #
## load("data/mc12artrejEEGLAB.rda")
## dat <- dat[dat$Time >= -200 & dat$Time <= 1500, ]
## rownames(dat) <- 1:nrow(dat)
## #
## # get blinks
## peaks.eeglab <- get.peaks(dat,"E21",NULL)
## save(peaks.eeglab, file = "data/peaks.eeglab.rda", compress= "xz")
## #
## # insert event code 777 at each peak
## dat$EventCode <- as.character(dat$EventCode)
## pb<-txtProgressBar(min=1,max=length(peaks.eeglab),char="=",
## 	style=3)
## for(i in 1:length(peaks.eeglab)){
## 	setTxtProgressBar(pb,i)
## 	tmp <- peaks.eeglab[[i]]
## 	if(!is.na(tmp[1])){
## 		for(j in 1:length(tmp)){
## 			dat[dat$Trial==i & dat$Time==tmp[j], 
## 				"EventCode"] <- "777"
## 		}
## 	}
## }
## close(pb)
## #
## # save event codes to later merge with corrected data frame
## evts <- dat$EventCode
## save(evts, file = "data/evts.peaks.eeglab.rda", compress = "xz")
## #
## # grab a 200 ms window around each peak and
## # put into data frame
## x <- as.numeric(rownames(dat[dat$EventCode=="777",]))
## x1 <- x - 25
## x2 <- x + 25
## x <- cbind(x1, x2)
## tmp <- dat[x[1, 1]:x[1, 2], ]
## tmp[1,]$EventCode <- "111111"
## tmp[nrow(tmp),]$EventCode <- "222222"
## pb <- txtProgressBar(min = 1, max = nrow(x), char = "=",
## 	style = 3)
## for(i in 2:nrow(x)){
## 	setTxtProgressBar(pb, i)
## 	tmp1 <- dat[x[i, 1]:x[i, 2], ]
## 	tmp1[1,]$EventCode <- "111111"
## 	tmp1[nrow(tmp1),]$EventCode <- "222222"
## 	tmp <- rbind(tmp, tmp1)
## }
## close(pb)
## rownames(tmp) <- 1:nrow(tmp)
## #
## table(tmp$EventCode)
## # For some reason, two 777 codes appear in some trials.
## # Remove the ones after the first 777.
## tmp<-add.trial2(tmp,list(begin=111111,finish=222222))
## pb <- txtProgressBar(min = 1, max = length(unique(tmp$Trial)), 
## 	char = "=", style = 3)
## for(i in 1:length(unique(tmp$Trial))){
## 	setTxtProgressBar(pb, i)
## 	while(length(which(tmp[tmp$Trial==i,"EventCode"]==777))>1){
## 		tmp[as.numeric(rownames(tmp[tmp$Trial==i&tmp$EventCode==777,
## 			]))[2],"EventCode"]<-0
## 	}
## }
## close(pb)
## #
## # reset time with t = 0 at event code 777
## tmp <- add.time2(x = tmp, markers = list(begin = "111111", 
##      ref = "777", finish = "222222"), sampling.rate = 250)
## eeglab.uncor <- tmp
## save(eeglab.uncor, file = "data/eeglab.uncor.rda", compress = "xz")
## rm(tmp); gc(TRUE, TRUE)
## #
## # compute blink average for corrected data
## datc <- eeglab$data
## datc$EventCode <- evts
## #
## # grab a 200 ms window around each peak and
## 
## # put into data frame
## x <- as.numeric(rownames(datc[datc$EventCode=="777",]))
## x1 <- x - 25
## x2 <- x + 25
## x <- cbind(x1, x2)
## tmp <- datc[x[1, 1]:x[1, 2], ]
## tmp[1,]$EventCode <- "111111"
## tmp[nrow(tmp),]$EventCode <- "222222"
## pb <- txtProgressBar(min = 1, max = nrow(x), char = "=",
## 	style = 3)
## for(i in 2:nrow(x)){
## 	setTxtProgressBar(pb, i)
## 	tmp1 <- datc[x[i, 1]:x[i, 2], ]
## 	tmp1[1,]$EventCode <- "111111"
## 	tmp1[nrow(tmp1),]$EventCode <- "222222"
## 	tmp <- rbind(tmp, tmp1)
## }
## rownames(tmp) <- 1:nrow(tmp)
## #
## table(tmp$EventCode)
## # For some reason, two 777 codes appear in some trials.
## # Remove the ones after the first 777.
## tmp<-add.trial2(tmp,list(begin=111111,finish=222222))
## pb <- txtProgressBar(min = 1, max = length(unique(tmp$Trial)), 
## 	char = "=", style = 3)
## for(i in 1:length(unique(tmp$Trial))){
## 	setTxtProgressBar(pb, i)
## 	while(length(which(tmp[tmp$Trial==i,"EventCode"]==777))>1){
## 		tmp[as.numeric(rownames(tmp[tmp$Trial==i&tmp$EventCode==777,
## 			]))[2],"EventCode"]<-0
## 	}
## }
## close(pb)
## #
## # reset time with t = 0 at event code 777
## tmp <- add.time2(x = tmp, markers = list(begin = "111111", 
##      ref = "777", finish = "222222"), sampling.rate = 250)
## eeglab.cor <- tmp
## save(eeglab.cor, file = "data/eeglab.cor.rda", compress = "xz")
## rm(tmp); gc(TRUE, TRUE)
## #
## # topomap for uncorrected
## avg <- tapply(eeglab.uncor[,1], eeglab.uncor$Time, mean)
## chan<-colnames(dat)[1:129]
## avg.dat<-data.frame(Time=as.numeric(names(avg)),Amplitude=avg,Channel=chan[1])
## pb <- txtProgressBar(min = 1, max = length(chan)-1, char = "=",
## 	style = 3)
## for(i in 2:length(chan)){
## 	setTxtProgressBar(pb, i)
## 	avg <- tapply(eeglab.uncor[,chan[i]], eeglab.uncor$Time, mean)
## 	avg.dat<-rbind(avg.dat,data.frame(Time=as.numeric(names(avg)),
## 		Amplitude=avg,Channel=chan[i]))
## }
## close(pb)
## avg.dat$Channel<-as.factor(avg.dat$Channel)
## coords<-des("egi.129")$cart
## avg.eeglab.uncor<-merge(avg.dat,coords[,c("x","y","Channel")],by="Channel")
## m.uncor<-gam(Amplitude~te(x,y,bs="ts",k=11),dat=avg.eeglab.uncor)
## #
## # topomap for corrected
## avg <- tapply(eeglab.cor[,1], eeglab.cor$Time, mean)
## chan<-colnames(dat)[1:129]
## avg.dat<-data.frame(Time=as.numeric(names(avg)),Amplitude=avg,Channel=chan[1])
## pb <- txtProgressBar(min = 1, max = length(chan)-1, char = "=",
## 	style = 3)
## for(i in 2:length(chan)){
## 	setTxtProgressBar(pb, i)
## 	avg <- tapply(eeglab.cor[,chan[i]], eeglab.cor$Time, mean)
## 	avg.dat<-rbind(avg.dat,data.frame(Time=as.numeric(names(avg)),
## 		Amplitude=avg,Channel=chan[i]))
## }
## close(pb)
## avg.dat$Channel<-as.factor(avg.dat$Channel)
## coords<-des("egi.129")$cart
## avg.eeglab.cor<-merge(avg.dat,coords[,c("x","y","Channel")],by="Channel")
## m.cor<-gam(Amplitude~te(x,y,bs="ts",k=11),dat=avg.eeglab.cor)
## #
## # get plotting info
## pi.eeglab.uncor<-plotGAM(m.uncor,too.far=des("egi.129")$too.far,plot=FALSE)
## pi.eeglab.cor<-plotGAM(m.cor,too.far=des("egi.129")$too.far,plot=FALSE)
## #
## # waveforms
## avg.uncor <- tapply(eeglab.uncor$E21, eeglab.uncor$Time, mean)
## avg.cor <- tapply(eeglab.cor$E21, eeglab.cor$Time, mean)
## time <- as.numeric(names(avg.uncor))
## #
## # create plot
## pdf(file = "EEGLABuncorrectedCorrected.pdf")
## par(mfrow = c(2, 2), mar = c(5.1, 4.1, 4.1, 2.1))
## # plot waveforms
## plot(time, avg.uncor - min(avg.uncor), type = "l",
## 	xlab = "Time (ms)", ylab = "Amplitude",ylim=c(0,500),
## 	main = "E21")
## lines(time, avg.cor - min(avg.cor), col = 2)
## legend("topleft", legend = c("blinks in uncorrected data",
## 	"blinks in corrected data"), lty = 1, col = 1:2, bty = "n")
## par(mar = c(1, 1, 1, 1))
## zlimit <- range(rbind(pi.eeglab.uncor$mat, pi.eeglab.cor$mat), na.rm = TRUE)
## # skip a plotting region
## plot.new()
## # plot uncorrected topomap
## image(pi.eeglab.uncor$xm, pi.eeglab.uncor$xm, pi.eeglab.uncor$mat, col = topo.colors(100),
## 	xlab = "", ylab = "", zlim = zlimit, axes = FALSE)
## contour(pi.eeglab.uncor$xm, pi.eeglab.uncor$xm, pi.eeglab.uncor$mat, add = TRUE)
## title(main = "Uncorrected", line = -12.5)
## # plot corrected topomap
## image(pi.eeglab.cor$xm, pi.eeglab.cor$xm, pi.eeglab.cor$mat, col = topo.colors(100),
## 	xlab = "", ylab = "", zlim = zlimit, axes = FALSE)
## contour(pi.eeglab.cor$xm, pi.eeglab.cor$xm, pi.eeglab.cor$mat, add = TRUE)
## title(main = "Corrected", line = -12.5)
## par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
## dev.off()
## #
## save(avg.eeglab.cor,avg.eeglab.uncor,eeglab.cor,eeglab.uncor,
## 	pi.eeglab.cor,pi.eeglab.uncor,m.uncor,m.cor,
## 	file="models/eeglab_correction.rda",compress="xz")


###################################################
### code chunk number 33: EEGLABvsICAocularCorrection (eval = FALSE)
###################################################
## pdf(file = "EEGLABvsICAocularCorrection.pdf")
## par(mfrow = c(2, 2), mar = c(5.1, 4.1, 4.1, 2.1))
## # plot waveforms
## avg.uncor <- tapply(eeglab.uncor$E21, eeglab.uncor$Time, mean)
## avg.cor <- tapply(eeglab.cor$E21, eeglab.cor$Time, mean)
## time <- as.numeric(names(avg.uncor))
## plot(time[1:51], avg.uncor[1:51] - min(avg.uncor[1:51]), type = "l",
## 	xlab = "Time (ms)", ylab = "Amplitude",ylim=c(0,600),
## 	main = "E21")
## lines(time[1:51], avg.cor[1:51] - min(avg.cor[1:51]), col = 2)
## avg.uncor <- tapply(eog.uncor$E21, eog.uncor$Time, mean)
## avg.cor <- tapply(eog.cor3$E21, eog.cor3$Time, mean)
## time <- as.numeric(names(avg.uncor))
## lines(time, avg.uncor - min(avg.uncor), col = 3)
## lines(time, avg.cor - min(avg.cor), col = 4)
## legend("topleft", legend = c("EEGLAB -- blinks in uncorrected data",
## 	"EEGLAB -- blinks in corrected data", 
## 	"icOC -- blinks in uncorrected data",
## 	"icOc -- blinks in corrected data"), 
## 	 lty = 1, col = 1:4, bty = "n",cex=0.85)
## par(mar = c(1, 1, 1, 1))
## zlimit <- range(rbind(pi.uncor$mat, pi.cor3$mat, 
## 	pi.eeglab.cor$mat), na.rm = TRUE)
## # uncorrected topomap
## image(pi.uncor$xm, pi.uncor$xm, pi.uncor$mat, col = topo.colors(100),
## 	xlab = "", ylab = "", zlim = zlimit, axes = FALSE)
## contour(pi.uncor$xm, pi.uncor$xm, pi.uncor$mat, add = TRUE)
## title(main = "icaOcularCorrection \n corrected", line = -12.5)
## #
## # icaOcularCorrection -- corrected topomap
## image(pi.cor3$xm, pi.cor3$xm, pi.cor3$mat, col = topo.colors(100),
## 	xlab = "", ylab = "", zlim = zlimit, axes = FALSE)
## contour(pi.cor3$xm, pi.cor3$xm, pi.cor3$mat, add = TRUE)
## title(main = "icaOcularCorrection \n corrected", line = -12.5)
## #
## # EEGLAB --  corrected topomap
## image(pi.eeglab.cor$xm, pi.eeglab.cor$xm, pi.eeglab.cor$mat, 
## 	col = topo.colors(100), xlab = "", ylab = "", zlim = zlimit, 
## 	axes = FALSE)
## contour(pi.eeglab.cor$xm, pi.eeglab.cor$xm, pi.eeglab.cor$mat, 
## 	add = TRUE)
## title(main = "EEGLAB \n Corrected", line = -12.5)
## par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
## dev.off()


