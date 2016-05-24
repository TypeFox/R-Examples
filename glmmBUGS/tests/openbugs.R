

if( "R2OpenBUGS" %in% rownames(installed.packages)) {
	
	
	
	library(MASS)
	data(bacteria)
	
	bacterianew <- bacteria
	bacterianew$yInt = as.integer(bacterianew$y == "y")
	levels(bacterianew$trt) <- c("placebo",
			"drug", "drugplus")
	
	library(glmmBUGS)
	bacrag <- glmmBUGS(formula = yInt ~ trt + week, 
			data = bacterianew,
			effects = "ID", modelFile = "bacteria.txt",
			family = "bernoulli",brugs=TRUE)
	
	
	source("getInits.R")
	startingValues = bacrag$startingValues
	
	library(R2OpenBUGS)
	
	bacResult = bugs(bacrag$ragged, inits=getInits,
			model.file = "bacteria.txt", n.chain = 3,
			n.iter = 200, n.burnin = 10,
			parameters = names(getInits()),
			n.thin = 4)
	
	bacParams = restoreParams(bacResult,
			bacrag$ragged)
	bacsummary = summaryChain(bacParams)
	
	
	bacsummary$betas
	
	
	checkChain(bacParams, c("intercept", "SDID"),oneFigure=FALSE)
	
	
	
	
	
}