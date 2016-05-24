IASD <-
function(df, dfCols = NA, fixSignApproximation = FALSE, 
		plotGraph = TRUE, plotToScreen = FALSE, filePrefix = NA, 
		xlimMin = NA, xlimMax = NA, ylimMin = 0, ylimMax = NA, 
		dHist = NA, dFunc = NA, meanStartSymmetric = NA, 
		sdStartSymmetric = NA, meanStartAsymmetric = NA, 
		sdStartAsymmetric = NA, positiveRatioStartAsymmetric = NA, 
		plotSelect = rep(TRUE, 4), showLegend = TRUE, 
		modelName = c("FA", "DA", "AS", "Skewed AS"), xlab = NA, 
		ylab = NA, main = NA, freqAxis = FALSE, lineColor = "black", 
		nsmall = 2, fileType = "TEXT", ...){  
	# df: data frame
	# range: range of the distribution graph, 
	# dHist: division width of histgram
	# dFunc: division of function graph
	
	#Sample Variance. var in R is unbiased variance
	BVar <- function(x){
		sum((x-mean(x))*(x-mean(x)))/length(x)
	}
	
	#used to calculated AICc, the small-sample-size corrected version of AIC
	AICC <- function(p, n){
		p*n/(n-p-1)
	}

	PlotHist <- function(ylimMin, ylimMax, freq = F){
		if (is.na(ylimMin) || is.na(ylimMax))
			hist(ias, freq = freq, breaks = breaks, xlab = xlab[i], ylab = ylab[i], main = main[i], ...)
		else
			hist(ias, freq = freq, breaks = breaks, ylim = c(ylimMin, ylimMax), xlab = xlab[i], ylab = ylab[i], main = main[i], ...)
	}
	
	AddPrefix <- function(name){
		if (is.na(filePrefix))
			return(name)
		else
			return(paste0(filePrefix, "-", name))
	}
		
	if (is.na(dfCols)){
		ncdf <- ncol(df)
		if (ncdf > 1)
			dfCols = c(2:ncdf)
		else
			dfCols <- 1
	}
	nc <- length(dfCols)
	sd1 <- 0
	sd2 <- 0
	sd3 <- 0
	sd4 <- 0
	mean2 <- 0
	mean3 <- 0
	mean4 <- 0
	positiveRatio4 <- 0
	UniSymAIC <- rep(0, nc)
	UniSymAICC <- rep(0, nc)
	UniAsymAIC <- rep(0, nc)
	UniAsymAICC <- rep(0, nc)
	BiSymAIC <- rep(0, nc)
	BiSymAICC <- rep(0, nc)
	BiAsymAIC <- rep(0, nc)
	BiAsymAICC <- rep(0, nc)
	UnimodalSymmetric <- as.list(NA, nc)
	UnimodalAsymmetric <- as.list(NA, nc)
	BimodalSymmetric <- as.list(NA, nc)
	BimodalAsymmetric <- as.list(NA, nc)
	
	AIC.df <- data.frame(dummy = rep(0,4))
	rownames(AIC.df) <- modelName
	AICc.df <- data.frame(dummy = rep(0,4))
	rownames(AICc.df) <- modelName
	
	RepVar <- function(x){
		if (length(x) >= nc)
			return(x)
		else
			return(rep(x, nc))
	}
	if (nc > 1){
		fixSignApproximation <- RepVar(fixSignApproximation)
		xlimMin <- RepVar(xlimMin)
		xlimMax <- RepVar(xlimMax)
		dHist <- RepVar(dHist)
		dFunc <- RepVar(dFunc)
		ylimMin <- RepVar(ylimMin)
		ylimMax <- RepVar(ylimMax)
		meanStartSymmetric <- RepVar(meanStartSymmetric)
		sdStartSymmetric <- RepVar(sdStartSymmetric)
		meanStartAsymmetric <- RepVar(meanStartAsymmetric)
		sdStartAsymmetric <- RepVar(sdStartAsymmetric)
		positiveRatioStartAsymmetric <- RepVar(positiveRatioStartAsymmetric)
		xlab <- RepVar(xlab)
		ylab <- RepVar(ylab)
		main <- RepVar(main)
	}
	if (length(lineColor) == 1)
		lineColor <- rep(lineColor, 4)
	i <- 0
	for(clm in dfCols){
		i <- i+1
		lty <- 0
		colname <- colnames(df)[clm]
		ias <- df[[clm]]
		n <- length(ias)
		absMax = max(abs(ias))
		if (absMax > 6)
			absMax <- ceiling(absMax/5)*5
		else
			absMax <- ceiling(absMax)
		if (is.na(xlimMin[i]))
			xlimMin[i] <- - absMax
		if (is.na(xlimMax[i]))
			xlimMax[i] <- absMax
		if (is.na(dHist[i]))
			dHist[i] <- (xlimMax[i] - xlimMin[i])/20
		if (is.na(dFunc[i]))
			dFunc[i] <- (xlimMax[i] - xlimMin[i])/200
		if (plotGraph){
			breaks <- seq(xlimMin[i], xlimMax[i], dHist[i])
			iv <- seq(xlimMin[i], xlimMax[i], dFunc[i])
		}
		sd <- sqrt(sum(ias*ias)/n)
		df$PUniSym <- dnorm(ias, sd = sd)
		UniSymAIC[i] <- -2*sum(log(df$PUniSym)) + 2*1
		UniSymAICC[i] <- -2*sum(log(df$PUniSym)) + 2*AICC(1,n)

		f <- eval(parse(text = paste("function(x){dnorm(x, sd =", sd, ")}", 
			srcfile = NULL)))	

		UnimodalSymmetric[[i]] <- list(AIC = UniSymAIC[i], AICc = UniSymAICC[i], 
			sd = sd, f = f)
		names(UnimodalSymmetric)[i] = colname
		if (plotGraph){
			if (!plotToScreen)
				pdf(file = AddPrefix(paste0(colname, ".pdf")))
			if (is.na(xlab[i]))
				xlab[i] <- colname
			if (is.na(main[i]))
				main[i] <- sprintf("Histogram of %s", colname)
			if (is.na(ylab[i]))
				ylab[i] <- "Density"
			hr <- PlotHist(ylimMin[i], ylimMax[i])
			if (plotSelect[1]){
				lty <- lty + 1
				if (length(lineColor) == 1 || lineColor[1] == lineColor[2])
					ltyl <- lty
				else
					ltyl <- 1
				lines(iv, dnorm(iv, sd = sd), lty = ltyl, col = lineColor[lty])
			}
		}
			
		mean <- mean(ias)
		sd <- sqrt(BVar(ias))
		df$PUniAsymB <- dnorm(ias, mean = mean, sd = sd)
		UniAsymAIC[i] <- -2*sum(log(df$PUniAsymB)) + 2*2
		UniAsymAICC[i] <- -2*sum(log(df$PUniAsymB)) + 2*AICC(2,n)	

		f <- eval(parse(text = paste("function(x){dnorm(x, mean =", mean, 
			", sd =", sd,")}", srcfile = NULL)))	

		UnimodalAsymmetric[[i]] = list(AIC = UniAsymAIC[i], AICc = UniAsymAICC[i], 
			mean = mean, sd = sd, f = f)
		names(UnimodalAsymmetric)[i] = colname
		if (plotGraph & plotSelect[2]){
			lty <- lty + 1
			if (length(lineColor) == 1 || lineColor[1] == lineColor[2])
				ltyl <- lty
			else
				ltyl <- 1
			lines(iv, dnorm(iv, mean = mean, sd = sd), lty = ltyl, 
				col = lineColor[lty])
		}
		
		DBiNorm <- function(ias, mean = 0, sd = 1, positiveRatio = 0.5) {
			positiveRatio*dnorm(ias, mean = mean, sd = sd) +
				 (1-positiveRatio)*dnorm(ias, mean = -mean, sd = sd)
		}
		mean <- mean(abs(ias))
		sd <- sqrt(BVar(abs(ias)))
		if (fixSignApproximation[i]){
			df$PBiB <- DBiNorm(ias, mean = mean, sd = sd)
			BiSymAIC[i] <- -2*sum(log(df$PBiB)) + 2*2
			BiSymAICC[i] <- -2*sum(log(df$PBiB)) + 2*AICC(2,n)		
		}else{
			if (is.na(meanStartSymmetric[i]))
				meanStartSymmetric[i] <- mean
			if (is.na(sdStartSymmetric[i]))
				sdStartSymmetric[i] <- sd
			BiSymFunc <- function(mean = 2, sd = 2)
					-sum(log(DBiNorm(ias, mean = mean, sd = sd)))			
			BiSymFit <- mle(BiSymFunc, start = list(mean = meanStartSymmetric[i], 
				sd = sdStartSymmetric[i]))						
			BiSymAIC[i] <- -2*as.numeric(logLik(BiSymFit)) + 2*2
			BiSymAICC[i] <- -2*as.numeric(logLik(BiSymFit)) + 2*AICC(2,n)			
			mean = coef(BiSymFit)["mean"]
			sd = coef(BiSymFit)["sd"]
		}

		f <- eval(parse(text = paste("function(x){0.5*dnorm(x, mean =", mean,
			", sd =", sd,") + 0.5*dnorm(x, mean =", - mean,", sd =", sd, 
			")}", srcfile = NULL)))	

		BimodalSymmetric[[i]] <- list(AIC = BiSymAIC[i], AICc = BiSymAICC[i], 
			mean = mean, sd = sd, f = f)
		names(BimodalSymmetric)[i] = colname
		if (plotGraph & plotSelect[3]){
			lty <- lty + 1
			if (length(lineColor) == 1 || lineColor[1] == lineColor[2])
				ltyl <- lty
			else
				ltyl <- 1
			lines(iv, DBiNorm(iv, mean = mean, sd = sd), lty = ltyl, 
				col = lineColor[lty])
		}	
		
		mean <- mean(abs(ias))
		sd <- sqrt(BVar(abs(ias)))
		positiveRatio <- sum(ias > 0)/n
		if (fixSignApproximation[i]){
			df$PBiAsymB <- DBiNorm(ias, mean = mean, sd = sd, positiveRatio = positiveRatio)
			BiAsymAIC[i] <- -2*sum(log(df$PBiAsymB)) + 2*3
			BiAsymAICC[i] <- -2*sum(log(df$PBiAsymB)) + 2*AICC(3,n)
		}else{
			if (is.na(meanStartAsymmetric[i]))
				meanStartAsymmetric[i] <- mean
			if (is.na(sdStartAsymmetric[i]))
				sdStartAsymmetric[i] <- sd
			if (is.na(positiveRatioStartAsymmetric[i]))
				positiveRatioStartAsymmetric[i] <- positiveRatio
			BiAsymFunc <- function(mean = 2, sd =1 , positiveRatio = 0.5){
				if(positiveRatio < 0)
					positiveRatio <- 0
				if(positiveRatio > 1)
					positiveRatio <- 1
				-sum(log(DBiNorm(ias, mean = mean, sd = sd, positiveRatio = positiveRatio)))
			}
			#BiAsymFit <- mle(BiAsymFunc,lower=c(-100,0,0), upper=c(100,50,1)) 
				#bounds can only be used with method L-BFGS-B or Brent
			BiAsymFit <- mle(BiAsymFunc, start = list(mean = meanStartAsymmetric[i], 
				sd = sdStartAsymmetric[i], positiveRatio = positiveRatioStartAsymmetric[i]))
			BiAsymAIC[i] <- -2*as.numeric(logLik(BiAsymFit)) + 2*3
			BiAsymAICC[i] <- -2*as.numeric(logLik(BiAsymFit)) + 2*AICC(3,n)
			mean <- coef(BiAsymFit)["mean"]
			sd <- coef(BiAsymFit)["sd"]
			positiveRatio <- coef(BiAsymFit)["positiveRatio"]
		}

		f <- eval(parse(text = paste("function(x){", positiveRatio, "*dnorm(x, mean =", mean,
			", sd =", sd,") +", 1 - positiveRatio, "*dnorm(x, mean =", - mean,", sd =", sd, 
			")}", srcfile = NULL)))	

		BimodalAsymmetric[[i]] <- list(AIC = BiAsymAIC[i], AICc = BiAsymAICC[i], 
			mean = mean, sd = sd, positiveRatio = positiveRatio, f = f)
		names(BimodalAsymmetric)[i] = colname
		if (plotGraph){
			if (plotSelect[4]){
				lty <- lty + 1
				if (length(lineColor) == 1 || lineColor[1] == lineColor[2])
					ltyl <- lty
				else
					ltyl <- 1
				lines(iv, DBiNorm(iv, mean = mean, sd = sd, 
						positiveRatio = positiveRatio), lty = ltyl, 
					col = lineColor[lty])
			}
			if (showLegend){
				if (length(lineColor) == 1 || lineColor[1] == lineColor[2])
					ltyl <- 1:sum(plotSelect)
				else
					ltyl <- 1
				legend("topright", lty = ltyl, legend = modelName[plotSelect], 
					col = lineColor)
			}
			if (freqAxis){
				par(new = T)
				if (is.na(ylimMax[i]))
					ylim = NULL
				else
					ylim = c(0, ylimMax[i]*length(ias)
						*(hr$breaks[2] - hr$breaks[1]))
				hist(ias, freq = T, breaks = breaks, 
					ylim = ylim, xaxt = "n", yaxt = "n", 
					xlab = "", ylab = "", main = "")
				axis(4)
				mtext("Number of individuals", side=4, line = 3)
			}
			if (plotToScreen)			
				dev.print(pdf, file = AddPrefix(paste0(colname, ".pdf")))
			else
				dev.off()
		}

		AIC.df[i] <- c(UniSymAIC[i], UniAsymAIC[i], BiSymAIC[i], BiAsymAIC[i])
		colnames(AIC.df)[i] <- colname
		AICc.df[i] <- c(UniSymAICC[i], UniAsymAICC[i], BiSymAICC[i], BiAsymAICC[i])
		colnames(AICc.df)[i] <- colname
	}
	
	cat("AIC:\n")
	print(format(AIC.df, digits = 2, nsmall = nsmall))
	cat("\n")
	cat("AICc:\n")
	print(format(AICc.df, digits = 2, nsmall = nsmall))
	cat("\n")
	
	if (fileType == "CSV"){
		write.csv(AIC.df, file = AddPrefix("AIC.csv"))
		write.csv(AICc.df, file = AddPrefix("AICc.csv"))
	}
	
	if (fileType == "TEXT"){
		colN <- colnames(AIC.df)[1]
		colnames(AIC.df)[1] <- paste0("\t", colN)
		write.table(format(AIC.df, digits = 2, nsmall = nsmall), col.names = T, 
			file = AddPrefix("AIC.txt"), append = F, quote = F, sep = "\t")
		colN <- colnames(AICc.df)[1]
		colnames(AICc.df)[1] <- paste0("\t", colN)
		write.table(format(AICc.df, digits = 2, nsmall = nsmall), col.names = T, 
			file = AddPrefix("AICc.txt"), append = F, quote = F, sep = "\t")
	}

	result <- list(AIC = AIC.df, AICc = AICc.df, UnimodalSymmetric = UnimodalSymmetric, 
		UnimodalAsymmetric = UnimodalAsymmetric,
		BimodalSymmetric = BimodalSymmetric, 
		BimodalAsymmetric = BimodalAsymmetric)
	names(result)[3:6] <- modelName
	return(result)
}
