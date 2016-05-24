# last modified 2010-02-24 by S. Weisberg

sigma.hat <- function (...) {
	.Deprecated("sigmaHat", package="alr3")
	sigmaHat(...)
}

conf.intervals <- function (...) {
	.Deprecated("confint", package="alr3")
	confint(...)
}

boot.case <- function (...) {
	.Deprecated("bootCase", package="alr3")
	car:::bootCase(...)
}

pure.error.anova <- function (...) {
	.Deprecated("pureErrorAnova", package="alr3")
	pureErrorAnova(...)
}

delta.method <- function (...) {
	.Deprecated("deltaMethod", package="alr3")
	car:::deltaMethod(...)
}

powtran <- function (...) {
	.Deprecated("bcPower", package="alr3")
	car:::bcPower(...)
}

inv.tran.plot <- function (...) {
	.Deprecated("invTranPlot", package="alr3")
	car:::invTranPlot(...)
}

inv.tran.estimate <- function (...) {
	.Deprecated("invTranEstimate", package="alr3")
	car:::invTranEstimate(...)
}

inverse.response.plot <- function (...) {
	.Deprecated("inverseResponsePlot", package="alr3")
        car:::invResPlot(...)
}

inv.res.plot <- function (...) {
	.Deprecated("invResPlot", package="alr3")
	car:::invResPlot(...)
}

bctrans <- function (...) {
	.Deprecated("powerTransform", package="alr3")
	car:::powerTransform(...)
}

bctrans1 <- function (...) {
	.Deprecated("powerTransform", package="alr3", 
	"'bctrans1' is deprecated.\nUse 'powerTransform' in the 'car' package.\nThe arguments may have changed.")
	car:::powerTransform(...)
}

lrt.bctrans <- function (...) {
	.Deprecated("testTransform", package="alr3")
	car:::testTransform(...)
}

resid.curv.test <- function (...) {
	.Deprecated("residCurvTest", package="alr3")
	car:::residCurvTest(...)
}

tukey.nonadd.test <- function (...) {
	.Deprecated("tukeyNonaddTest", package="alr3")
	car:::tukeyNonaddTest(...)
}

resplot <- function (...) {
	.Deprecated("residualPlot", package="alr3")
	car:::residualPlot(...)
}

residual.plots <- function (...) {
	.Deprecated("residualPlots", package="alr3")
	car:::residualPlots(...)
}

marginal.model.plot <- function (...) {
	.Deprecated("marginalModelPlot", package="alr3")
	car:::marginalModelPlot(...)
}

marginal.model.plots <- function (...) {
	.Deprecated("marginalModelPlots", package="alr3")
	car:::marginalModelPlot(...)
}

inf.index <- function (...) {
	.Deprecated("infIndexPlot", package="alr3")
	car:::infIndexPlot(...)
}

outlier.t.test <- function (...) {
	.Deprecated("outlierTest", package="alr3")
	car:::outlierTest(...)
}

