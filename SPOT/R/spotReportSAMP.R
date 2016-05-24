###################################################################################################
#' SAMP Model Report
#' 
#' This function does a random effects model analysis of a SPOT run, using the lme4 package. It is supposed to be used
#' for the case of Single Algorithm Multiple Problems. The SPOT demo 21 provides a simple example for mixed model analysis.
#' Call for the demo is: \code{demo(spotDemo21MixedModelSAMP,ask=F)}.
#'
#' @param spotConfig the configuration list of all spot parameters
#' @return list spotConfig with changed values
#' @seealso  \code{\link{spotReportMAMP}}  \code{\link{spotStepReport}} 
#' @export
###################################################################################################
spotReportSAMP <- function(spotConfig) {	
	spotWriteLines(spotConfig$io.verbosity,2,"  Entering spotReportSAMP")
	
	######
	spotInstAndLoadPackages(c("lme4","ggplot2"))	#TODO in suggest		
	rawB <- spotConfig$alg.currentResult

	if(any(rawB$Y < 0))	
		yLog= log(rawB$Y-min(rawB$Y)+1)
	else
		yLog= log(rawB$Y)
	samp.df <- data.frame(cbind(y=rawB$Y, yLog=yLog, algSeed = rawB$SEED, fSeed=rawB$PINST))
	samp.df$algSeed <- factor(samp.df$algSeed)
	samp.df$fSeed <- factor(samp.df$fSeed)
	# 	samp.aov <- aov(yLog ~ fSeed, data=samp.df)
	# 	(M1 <- anova(samp.aov))
	# 	(MSA <- M1[1,3])
	# 	(MSE <- M1[2,3])
	# 	r <-length(unique(samp.df$algSeed))
	# 	q <- nlevels(samp.df$fSeed)
	#   print("1st variance component (treatment):")
	# 	(var.A <- (MSA - MSE)/(r))
	#   print("2nd variance component (error):")
	# 	(var.E <- MSE)
	# 	var.A + var.E
	# 	coef(samp.aov)[1]
	# 	1-pf(MSA/MSE,q-1,q*(r-1))
	# 	MSA.anova <- MSA
	###	
	samp.lmer <- lme4::lmer(y~ 1 +(1|fSeed),data=samp.df)
	spotPrint(spotConfig$io.verbosity,1,paste("Summary of the mixed model: ",sep=""));
	spotPrint(spotConfig$io.verbosity,1,samp.lmer)
	
  ### NEW: Check Model Adequacy:
	# checking that e_i all have the same variance:
	#pdf(file="qq1.pdf")
	dev.new()
	plot(resid(samp.lmer) ~ fitted(samp.lmer),main="residual plot")
	abline(h=0)
	#dev.off()
  
	dev.new()
	# checking the normality of residuals e_ij:
	qqnorm(resid(samp.lmer), main="Q-Q plot for residuals")
	qqline(resid(samp.lmer))
	
  ### 
  ### Since the raw data model is are not adequate, we perform a log transformation: 	
	samp.lmer.log <- lme4::lmer(yLog~ 1 +(1|fSeed),data=samp.df)
	spotPrint(spotConfig$io.verbosity,1,paste("Summary of the mixed model (logY): ",sep=""));
	spotPrint(spotConfig$io.verbosity,1,samp.lmer.log)
	
	# checking that e_ij all have the same variance:
	#pdf(file="qq1Log.pdf")
	dev.new()
	plot(resid(samp.lmer.log) ~ fitted(samp.lmer.log),main="residual plot (log.)")
	abline(h=0)
	#dev.off()
	dev.new()
	# checking the normality of residuals e_ij:
	qqnorm(resid(samp.lmer.log), main="Q-Q plot for residuals (log.)")
	qqline(resid(samp.lmer.log))
	
  ### :WEN
	####
	VC <- lme4::VarCorr(samp.lmer.log)
	sigma.tau <- as.numeric(attr(VC$fSeed,"stddev"))
	sigma <- as.numeric(attr(VC,"sc"))
	q <- nlevels(samp.df$fSeed)
	r <- length(unique(samp.df$algSeed))
	MSA <- sigma^2+r*sigma.tau^2
	MSE <- sigma^2
	### diesen wert als p-wert ausgeben:
	pvalue=1-pf(MSA/MSE,q-1,q*(r-1))
	spotPrint(spotConfig$io.verbosity,1,paste("P-value log.: ",pvalue,sep=""));
	###	
	s <- sqrt(MSA/(q*r))
	Ydotdot <- mean(samp.df$yLog)
	#qsr <- qt(1-0.025,r)
	qsr <- qt(1-0.025,q*(r-1)) #Fix TBB, 30.04.2013
	### conf intervall ausgeben:
	#if(any(rawB$Y < 0))	
	#	confInt=c( exp(Ydotdot - qsr * s)+(min(rawB$Y)-1), exp(Ydotdot + qsr * s)+(min(rawB$Y)-1))
	#else
		confInt.log=c( (Ydotdot - qsr * s), (Ydotdot + qsr * s))
	
	##############################
	VC <- lme4::VarCorr(samp.lmer)
	sigma.tau <- as.numeric(attr(VC$fSeed,"stddev"))
	sigma <- as.numeric(attr(VC,"sc"))
	q <- nlevels(samp.df$fSeed)
	r <- length(unique(samp.df$algSeed))
	MSA <- sigma^2+r*sigma.tau^2
	MSE <- sigma^2
	### diesen wert als p-wert ausgeben:
	pvalue=1-pf(MSA/MSE,q-1,q*(r-1))
	spotPrint(spotConfig$io.verbosity,1,paste("P-value: ",pvalue,sep=""));
	###	
	s <- sqrt(MSA/(q*r))
	Ydotdot <- mean(samp.df$y)
	#qsr <- qt(1-0.025,r)
	qsr <- qt(1-0.025,q*(r-1)) #Fix TBB, 30.04.2013
		### conf intervall ausgeben:
	#if(any(rawB$Y < 0))	
	#	confInt=c( exp(Ydotdot - qsr * s)+(min(rawB$Y)-1), exp(Ydotdot + qsr * s)+(min(rawB$Y)-1))
	#else
		confInt=c( (Ydotdot - qsr * s), (Ydotdot + qsr * s))
	##############################
	
	
	

	spotPrint(spotConfig$io.verbosity,1,paste("Confidence Interval log.: ",confInt.log[1]," to ", confInt.log[2], sep=""));
	spotPrint(spotConfig$io.verbosity,1,paste("Confidence Interval: ",confInt[1]," to ", confInt[2], sep=""));
	### same analysis based on anova
	# 	s <- sqrt(MSA.anova/(q*r))
	# 	Ydotdot <- mean(samp.df$yLog)
	# 	qsr <- qt(1-0.025,r)
	# 	c( exp(Ydotdot - qsr * s), exp(Ydotdot + qsr * s))
	#### plots  
	# 	ggplot(samp.df, aes(x = y, y = fSeed, colour = algSeed)) +
	# 	  geom_point() + opts(title = "Performance")
	#   pdf(file="gg1Log.pdf")  
	### nebeneinander plotten:
	dev.new()
	plt1 <- ggplot2::ggplot(samp.df, ggplot2::aes_string(x = "yLog", y = "fSeed")) + ggplot2::geom_point() + ggplot2::ggtitle("Performance")
	print(plt1)
	dev.new()
	plt2 <- ggplot2::ggplot(samp.df, ggplot2::aes_string(x = "y", y = "fSeed")) + ggplot2::geom_point() + ggplot2::ggtitle("Performance")
	print(plt2)
	#dev.off()
	spotConfig
}
