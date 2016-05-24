
## ----echo=FALSE,message=FALSE,results='hide'-----------------------------
options(markdown.HTML.stylesheet = 'extra/manual.css')
library(knitr)
options(digits=3)
require(graphics)
set.seed(2)


## ----eval=FALSE----------------------------------------------------------
## install.packages("BayesSingleSub")


## ------------------------------------------------------------------------
library(BayesSingleSub)


## ------------------------------------------------------------------------
data = c(87.5, 82.5, 53.4, 72.3, 94.2, 96.6, 57.4, 78.1, 
         47.2, 80.7, 82.1, 73.7, 49.3, 79.3, 73.3, 57.3, 31.7, 
         50.4, 77.8, 67, 40.5, 1.6, 38.6, 3.2, 24.1)

n1 = 10
n2 = 15


## ------------------------------------------------------------------------
ypre =  data[1:n1]
ypost = data[n1 + 1:n2]


## ----results='hide'------------------------------------------------------
outMeanDiff = ttest.Gibbs.AR(ypre, ypost, 
    iterations = 10000, return.chains = TRUE,
		r.scale = 1, betaTheta = 5, 
    areaNull=c(-.2,.2), return.onesided = TRUE, leftSided = TRUE,
    sdMet = 0.3)
outTrendDiff = trendtest.Gibbs.AR(ypre, ypost,
		iterations = 10000, return.chains = TRUE,
		r.scaleInt = 1, r.scaleSlp = 1, betaTheta = 5, 
    intArea=c(-.2,.2), slpArea=c(-.2,.2), 
    return.onesided = TRUE, leftSidedInt = TRUE, leftSidedSlp = TRUE,
    sdMet = 0.3)


## ------------------------------------------------------------------------
logbfMeanDiff = outMeanDiff$logbf
logbfMeanDiff.area =  outMeanDiff$logbfArea
logbfMeanDiff.onesided = outMeanDiff$logbfOnesided
chainsMeanDiff = outMeanDiff$chains

logbfTrendDiff = outTrendDiff$logbf
logbfTrendDiff.area =  outTrendDiff$logbfArea
logbfTrendDiff.onesided = outTrendDiff$logbfOnesided
chainsTrendDiff = outTrendDiff$chains


## ------------------------------------------------------------------------
exp(logbfMeanDiff)
exp(logbfMeanDiff.area)
exp(logbfMeanDiff.onesided)

exp(logbfTrendDiff)
exp(logbfTrendDiff.area)
exp(logbfTrendDiff.onesided)


## ----eval=FALSE----------------------------------------------------------
## hist(chainsTrendDiff[, 8], main = "Posterior for autocorrelation coeff.",
##     xlim = c(0, 1))


## ----echo=FALSE----------------------------------------------------------
hist(chainsTrendDiff[, 8],main="Posterior for autocorrelation coeff.", xlim=c(0,1),xlab=expression(rho))


## ------------------------------------------------------------------------
summary(chainsTrendDiff[, 8])


## ----results='hide'------------------------------------------------------
logBAR = ttest.MCGQ.AR(ypre, ypost, 
    iterations = 10000,  r.scale = 1,  betaTheta = 5)
logBTRENDS = trendtest.MC.AR(ypre, ypost, 
    iterations = 10000, r.scaleInt = 1, r.scaleSlp = 1,  betaTheta = 5)


## ------------------------------------------------------------------------
logBAR
logBTRENDS
exp(logBAR)  
exp(logBTRENDS)


