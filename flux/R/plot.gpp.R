plot.gpp <- function(x, nm = "", single.pane = TRUE, ...){
	# prepare plotting region
	if(single.pane){par(mfrow=c(1,3), las=1, pty="s")}
	### plot only Reco
	plot(R ~ Temp, x$mr$model, main=paste(nm, "Reco model"), ...)
	t.range <- range(x$mr$model$Temp)
	ndr <- seq(min(t.range), max(t.range), length.out=100)
	lines(ndr, predict(x$mr, newdata=data.frame(Temp=ndr)))
	
	### GPP + NEE
	# first determine complete range
	flux.range <- range(x$mg$model$GPP, x$data$Reco, x$mg$model$GPP + x$data$dat$Reco + x$data$offset)*1.1
	plot(GPP ~ PAR, data=x$mg$model, ylab="GPP/NEE", ylim=flux.range, main = "Modeled NEE, GPP, Reco", ...)
	# NEE
	points(GPP + x$data$dat$Reco + x$data$offset ~ PAR, data=x$mg$model, pch="+", col="grey")
	# straight prediction
	PAR.ord <- order(x$mg$model$PAR)
	lines((predict(x$mg) + x$data$dat$Reco + x$data$offset)[PAR.ord] ~ x$mg$model$PAR[PAR.ord])
	# smooth prediction, therefore:
	# model of temperature vs PAR
	mpt <- loess(x$data$PAR.Temp ~ PAR, x$mg$model)
	# PAR new data for smooth lines
	nd <- with(x$mg$model, seq(min(PAR), max(PAR), length.out=100))
	# plot resultant NEE model
	lines(nd, predict(x$mg, newdata=data.frame(PAR=nd)) + predict(x$mr, newdata=data.frame(Temp = predict(mpt, newdata=data.frame(PAR=nd)))) + x$data$offset, col="green4")
	# plot GPP model alone
	lines(nd, predict(x$mg, newdata=data.frame(PAR=nd)), col="green4", lty=2)
	# plot Reco model alone
	lines(nd, predict(x$mr, newdata=data.frame(Temp = predict(mpt, newdata=data.frame(PAR=nd)))), col="green4", lty=3)
	# legend
	legend("bottomleft", bty="n", lty=c(1,2,3), col="green4", legend=c("NEE", "GPP", "Reco"))
	
	### NEE modelled vs measured with Fit
	nee.pred <- predict(x$mg) + x$data$dat$Reco + x$data$offset
	nee.meas <- x$mg$model$GPP + x$data$dat$Reco + x$data$offset
	nee.range <- range(nee.pred, nee.meas)
	plot(nee.meas ~ nee.pred, xlim=nee.range, ylim=nee.range, main="NEE measured vs predicted", ylab = "NEE measured", xlab = "NEE predicted", ...)
	abline(lm(nee.meas ~ nee.pred))
	abline(coef=c(0,1), lty=3)
	legend("topleft", bty="n", legend = c(paste("R2 = ", round(summary(lm(nee.meas ~ nee.pred))$r.squared, 2)), paste("Intercept = ", round(coef(lm(nee.meas ~ nee.pred))[1],3)), paste("Slope = ", round(coef(lm(nee.meas ~ nee.pred))[2],3))), cex=1.2)
}