# Test varia

require(MuMIn)
packageVersion("MuMIn")
options(na.action = "na.fail")

#print(packageDescription("MuMIn", fields = "Version"))
# TEST binary response ---------------------------------------------------------
ldose <- rep(0:5, 2)
numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
sex <- factor(rep(c("M", "F"), c(6, 6)))
SF <- cbind(numdead, numalive=20-numdead)
budworm.lg <- glm(SF ~ sex*ldose, family = binomial)
dd <- dredge(budworm.lg, trace=FALSE)
gm <- get.models(dd, 1:4)
model.avg(gm)

# The same, but use cbind directly in the formula
budworm.lg <- glm(cbind(numdead, numalive=20-numdead) ~ sex*ldose, family=binomial)
dd <- dredge(budworm.lg, trace=TRUE)
avgmod <- model.avg(get.models(dd, 1:4))

# TEST for consistency of vcov and se calculation ------------------------------
if(!isTRUE(all.equal(coefTable(avgmod, adjust.se = FALSE)[,2],
	sqrt(diag(vcov(avgmod))), tolerance = .001)))
	stop("'vcov' has a problem")

# TEST evaluation from within function -----------------------------------------

budworm <- data.frame(ldose = rep(0:5, 2), numdead = c(1, 4, 9, 12, 18, 20, 0,
	2, 6, 10, 12, 16), sex = factor(rep(c("M", "F"), c(6, 6))))
budworm$SF <- cbind(numdead = budworm$numdead, numalive = 20 - budworm$numdead)

# evaluate within an exotic environment
(function(dat) (function(dat2) {
	#mod <- glm(SF ~ sex*ldose, data = dat2, family = "quasibinomial", trace=T)
	mod <- glm(SF ~ sex*ldose, data = dat2, family = "binomial")
	#mod <- glm(SF ~ sex*ldose, data = budworm, family = "binomial", trace=F)
	print(dd <- dredge(mod, rank = "QAIC", chat = summary(budworm.lg)$dispersion))
	gm <- get.models(dd, subset = NA, family = "binomial")
	#print(sys.frames())
	summary(model.avg(gm))
})(dat))(budworm)



rm(list=ls())

# END TESTS
