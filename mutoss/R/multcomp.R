# Simultaneous confidence intervals for linear contrasts
# 
# Author: FrankKonietschke
###############################################################################







multcomp.wrapper <- function(model, hypotheses, alternative, rhs=0, alpha, factorC ) {
	
	type<-""
	
	if (any(factorC== c("Tukey", "Dunnett", "Sequen",
			"Williams", "Changepoint", "AVE", "McDermott", "Marcus", "UmbrellaWilliams"))) {
		eval(parse(text=paste("type <- mcp(",factorC,"=hypotheses)")))
	} else {		
		eval(parse(text=paste("type <- mcp(",factorC,"=hypotheses)")))
	}	

	glhtObj <- glht(model, linfct = type, rhs=rhs, alternative=alternative)
	summaryGLHT <- summary(glhtObj)
	pvalues <- summaryGLHT$test$pvalues
	estimates <- summaryGLHT$test$coefficients
	confint <- confint(glhtObj,level=(1-alpha))$confint
	
	
	
	rejected1 <- (pvalues < alpha)
	confi <- cbind(confint)
	print(cbind(confi,pvalues))
	return(list(adjPValues=pvalues,rejected=rejected1,confIntervals= confi,
					errorControl = new(Class='ErrorControl',type="FWER",alpha=alpha)))

	

}


mutoss.multcomp<- function() { return(new(Class="MutossMethod",
					label="Multiple Contrast Tests",
					errorControl="FWER",
					callFunction="multcomp.wrapper",
					output=c("adjPValues", "rejected","confIntervals","errorControl"),
					info="<h2>Parametric multiple contrast tests and simultaneous confidence intervals</h2>
							<p>With this function, it is possible to compute simultaneous tests and confidence intervals for general linear hypotheses in parametric models<p> 
							<p></p>
							<h3>Reference:</h3>
							<ul>
							<li>Frank Bretz, Alan Genz and Ludwig A. Hothorn \"<i>On the numerical availability of multiple comparison procedures.</i>\" Biometrical Journal, 43(5), 645-656. , 2001.</li>
													</ul>",
					parameters=list(model=list(type="ANY"),
							hypotheses=list(type="ANY"),
							alpha=list(type="numeric"),
							alternative=list(type="character", label="Alternative", choices=c("two.sided", "less", "greater")),
							factorC=list(type="character", label="Factor for Comparison", fromR="FactorVar")
						
					)
			)) }






