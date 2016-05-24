user.prompt <- function (x) {
 	ANSWER <- readline(paste("\nType 'y' to continue to", x,"or 'n' to quit: ",sep=" "))  
	if (substr(ANSWER, 1, 1) == "n")
         {stop("Function terminated by user.")}
	}

LupPon <- function()
{
    message("Replication of Table 2 of Lupu and Pontusson (2011)'...")
	
	cat("\n")
    cat("Loading data...\n")
    data(LupPon)
    cat("> data(LupPon)\n")
    
    cat("Subsetting data and calculating lagged DV...\n")
	cat("> LupPon <- LupPon[!is.na(LupPon$redist),]\n> LupPon$redist.lag <- unlist(by(LupPon,LupPon$id,function(x){c(NA,x[,\"redist\"][1:(length(x[,\"redist\"])-1)])}))\n> LupPon$time <- unlist(by(LupPon,LupPon$id,function(x) seq(1:nrow(x))))")
	LupPon <- LupPon[!is.na(LupPon$redist),]
	LupPon$redist.lag <- unlist(by(LupPon,LupPon$id,function(x){c(NA,x[,"redist"][1:(length(x[,"redist"])-1)])}))
	LupPon$time <- unlist(by(LupPon,LupPon$id,function(x) seq(1:nrow(x))))

	# Specification 1 
	user.prompt("specification 1")
	cat("> out1 <- panelAR(redist ~ redist.lag + ratio9050 + ratio5010 + turnout + fempar + propind + pvoc + union + unempl, data=LupPon, panelVar='id', timeVar='time', autoCorr='ar1', panelCorrMethod='pcse',rho.na.rm=TRUE, panel.weight='t-1', bound.rho=TRUE)\n")
	out1 <- panelAR(redist ~ redist.lag + ratio9050 + ratio5010 + turnout + fempar + propind + pvoc + union + unempl, data=LupPon, panelVar='id', timeVar='time', autoCorr='ar1', panelCorrMethod='pcse',rho.na.rm=TRUE, panel.weight='t-1', bound.rho=TRUE)
	cat("> summary(out1)\n")
	print(summary(out1))
	
	# Specification 2 (remove outliers)
	user.prompt("specification 2")
	cat("Removing outliers...\n")
	cat("> mod1.resid <- out1$residuals\nindex <- which(abs((mod1.resid-mean(mod1.resid))/sd(mod1.resid)) <= 1.5)\n> LupPon.nooutlier <- out1$model[index,]")
	mod1.resid <- out1$residuals
	index <- which(abs((mod1.resid-mean(mod1.resid))/sd(mod1.resid)) <= 1.5)
	LupPon.nooutlier <- out1$model[index,]
	
	cat("> out2 <- panelAR(redist ~ redist.lag + ratio9050 + ratio5010 + turnout + fempar + propind + pvoc + union + unempl, data=LupPon.nooutlier, panelVar='id', timeVar='time', autoCorr='ar1', panelCorrMethod='pcse',rho.na.rm=TRUE, panel.weight='t-1', bound.rho=TRUE)\n")
	out2 <- panelAR(redist ~ redist.lag + ratio9050 + ratio5010 + turnout + fempar + propind + pvoc + union + unempl, data=LupPon.nooutlier, panelVar='id', timeVar='time', autoCorr='ar1', panelCorrMethod='pcse',rho.na.rm=TRUE, panel.weight='t-1', bound.rho=TRUE)
	cat("> summary(out2)\n")
	print(summary(out2))
	
	# Specification 3
	user.prompt("specification 3")
	cat("> out3 <- panelAR(redist ~ ratio9050 + ratio5010 + as.factor(id), data=LupPon, panelVar='id', timeVar='time', autoCorr='ar1', panelCorrMethod='pcse',rho.na.rm=TRUE, panel.weight='t-1', bound.rho=TRUE)\n")
	out3 <- panelAR(redist ~ ratio9050 + ratio5010 + as.factor(id), data=LupPon, panelVar='id', timeVar='time', autoCorr='ar1', panelCorrMethod='pcse',rho.na.rm=TRUE, panel.weight='t-1', bound.rho=TRUE)
	cat("> summary(out3)\n")
	print(summary(out3))

	# Specification 4
	user.prompt("specification 4")
	cat("Removing outliers...\n")
	cat("> mod3.resid <- out3$residuals\nindex <- which(abs((mod3.resid-mean(mod3.resid))/sd(mod3.resid)) <= 1.5)\n> LupPon.nooutlier <- out3$model[index,]")
	mod3.resid <- out3$residuals
	index <- which(abs((mod3.resid-mean(mod3.resid))/sd(mod3.resid)) <= 1.5)
	LupPon.nooutlier <- out3$model[index,]
	cat("> out4 <- panelAR(redist ~ ratio9050 + ratio5010 + as.factor(id), data=LupPon.nooutlier, panelVar='id', timeVar='time', autoCorr='ar1', panelCorrMethod='pcse',rho.na.rm=TRUE, panel.weight='t-1', bound.rho=TRUE)\n")
	out4 <- panelAR(redist ~ ratio9050 + ratio5010 + as.factor(id), data=LupPon.nooutlier, panelVar='id', timeVar='time', autoCorr='ar1', panelCorrMethod='pcse',rho.na.rm=TRUE, panel.weight='t-1', bound.rho=TRUE)
	cat("> summary(out4)\n")
	print(summary(out4))
	

	# Specification 5
	user.prompt("specification 5")
	cat("> out5 <- panelAR(redist ~ redist.lag + ratio9010 + skew + turnout + fempar + propind + pvoc + union + unempl, data=LupPon, panelVar='id', timeVar='time', autoCorr='ar1', panelCorrMethod='pcse',rho.na.rm=TRUE, panel.weight='t-1', bound.rho=TRUE)\n")
	out5 <- panelAR(redist ~ redist.lag + ratio9010 + skew + turnout + fempar + propind + pvoc + union + unempl, data=LupPon, panelVar='id', timeVar='time', autoCorr='ar1', panelCorrMethod='pcse',rho.na.rm=TRUE, panel.weight='t-1', bound.rho=TRUE)
	cat("> summary(out5)\n")
	print(summary(out5))

	# Specification 6
	user.prompt("specification 6")
	cat("Removing outliers...\n")
	cat("> mod5.resid <- out5$residuals\nindex <- which(abs((mod5.resid-mean(mod5.resid))/sd(mod5.resid)) <= 1.5)\n> LupPon.nooutlier <- out5$model[index,]")
	mod5.resid <- out5$residuals
	index <- which(abs((mod5.resid-mean(mod5.resid))/sd(mod5.resid)) <= 1.5)
	LupPon.nooutlier <- out5$model[index,]
	cat("> out6 <- panelAR(redist ~ redist.lag + ratio9010 + skew + turnout + fempar + propind + pvoc + union + unempl, data=LupPon.nooutlier, panelVar='id', timeVar='time', autoCorr='ar1', panelCorrMethod='pcse',rho.na.rm=TRUE, panel.weight='t-1', bound.rho=TRUE)\n")
	out6 <- panelAR(redist ~ redist.lag + ratio9010 + skew + turnout + fempar + propind + pvoc + union + unempl, data=LupPon.nooutlier, panelVar='id', timeVar='time', autoCorr='ar1', panelCorrMethod='pcse',rho.na.rm=TRUE, panel.weight='t-1', bound.rho=TRUE)
	cat("> summary(out6)\n")
	print(summary(out6))
	

	# Specification 7
	user.prompt("specification 7")
	cat("> out7 <- panelAR(redist ~ ratio9010 + skew + as.factor(id), data=LupPon, panelVar='id', timeVar='time', autoCorr='ar1', panelCorrMethod='pcse',rho.na.rm=TRUE, panel.weight='t-1', bound.rho=TRUE)\n")
	out7 <- panelAR(redist ~ ratio9010 + skew + as.factor(id), data=LupPon, panelVar='id', timeVar='time', autoCorr='ar1', panelCorrMethod='pcse',rho.na.rm=TRUE, panel.weight='t-1', bound.rho=TRUE)
	cat("> summary(out7)\n")
	print(summary(out7))	

	# Specification 8
	user.prompt("specification 8")
	cat("Removing outliers...\n")
	cat("> mod7.resid <- out7$residuals\nindex <- which(abs((mod7.resid-mean(mod7.resid))/sd(mod7.resid)) <= 1.5)\n> LupPon.nooutlier <- out7$model[index,]")
	mod7.resid <- out7$residuals
	index <- which(abs((mod7.resid-mean(mod7.resid))/sd(mod7.resid)) <= 1.5)
	LupPon.nooutlier <- out7$model[index,]
	cat("> out8 <- panelAR(redist ~ ratio9010 + skew + as.factor(id), data=LupPon.nooutlier, panelVar='id', timeVar='time', autoCorr='ar1', panelCorrMethod='pcse',rho.na.rm=TRUE, panel.weight='t-1', bound.rho=TRUE)\n")
	out8 <- panelAR(redist ~ ratio9010 + skew + as.factor(id), data=LupPon.nooutlier, panelVar='id', timeVar='time', autoCorr='ar1', panelCorrMethod='pcse',rho.na.rm=TRUE, panel.weight='t-1', bound.rho=TRUE)
	cat("> summary(out8)\n")
	print(summary(out8))
}

LupPon()