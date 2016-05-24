
############## 
######## Examples
##############

x1 <- make.mtkFactor(name="x1", distribName="unif",
	 distribPara=list(min=-pi, max=pi))
x2 <- make.mtkFactor(name="x2", distribName="unif",
     distribPara=list(min=-pi, max=pi))
x3 <- make.mtkFactor(name="x3", distribName="unif", 
     distribPara=list(min=-pi, max=pi))
ishi.factors <- mtkExpFactors(list(x1,x2,x3))


Designer <- mtkNativeDesigner ("BasicMonteCarlo", 
		information=list(size=20))
model <- mtkNativeEvaluator("Ishigami" )
Analyser <- mtkNativeAnalyser("Regression", information=list(nboot=20) )

ishiReg <- new("mtkExpWorkflow", expFactors=ishi.factors,
		   processesVector=c(
				              design=Designer,
				              evaluate=model,
				              analyze=Analyser)
			  				)
run(ishiReg)
summary(ishiReg)
