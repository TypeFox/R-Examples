# Tukey HSD Test in parametric factorial designs
# 
# Author: FrankKonietschke
###############################################################################












tukey.wrapper <- function(model, alpha,factorC) {
	
	model <- aov(model)
	tukeyObj <- TukeyHSD(model,conf.level=1-alpha, factorC)
	

	estimates <- c(tukeyObj)[[1]][,1]
	confintL <- c(tukeyObj)[[1]][,2]
	confintU <- c(tukeyObj)[[1]][,3]
	confi <- cbind(estimates, confintL, confintU)
	pvalues <- c(tukeyObj)[[1]][,4]
	
	
	
	rejected1 <- (pvalues < alpha)
	#confi <- cbind(confint)
	print(cbind(confi,pvalues))
	return(list(adjPValues=pvalues,rejected=rejected1,confIntervals= confi,
					errorControl = new(Class='ErrorControl',type="FWER",alpha=alpha)))
}


#tukey.wrapper(aov(breaks~wool*tension,data=warpbreaks),alpha=0.05,"tension")

mutoss.tukey<- function() { return(new(Class="MutossMethod",
					label="Tukey HSD Test",
					errorControl="FWER",
					callFunction="tukey.wrapper",
					output=c("adjPValues", "rejected","confIntervals","errorControl"),
					info="<h2>Tukey HSD test and simultaneous confidence intervals in parametric factorial designs</h2>
							<p>'With this function, it is possible to compute all pairs comparisons for expectations and 
                              simultaneous confidence intervals in factorial linear models. Hereby, the all-pairs comparisons
                              can be performed for user given effects. The overall variance is estimated by the linear model
                              as well as the degree of freedom used by the studentized range distribution. <p> 
							<p></p>
							<h3>Reference:</h3>
							<ul>
							<li>Tukey, J.W. \"<i>The problem of multiple comparisons. Unpublished manuscript
reprinted in: The Collected Works of John W. Tukey (1953).</i>\" Volume 8 (1994),Braun, H.I. (ed.), Chapman and Hall, New York.</li>
							</ul>",
					parameters=list(model=list(type="ANY"),
							
							alpha=list(type="numeric"),
							factorC=list(type="character", label="Factor for Comparison", fromR="FactorVar")
					
					)
			)) }
