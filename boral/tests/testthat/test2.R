context("test auxilaryfunctions.R")


###############
true.lv <- rbind(rmvnorm(n=15,mean=c(1,2)),rmvnorm(n=15,mean=c(-3,-1))) 
lv.coefs <- cbind(matrix(runif(30*3),30,3),1)
X <- matrix(rnorm(30*4),30,4) 
X.coefs <- matrix(rnorm(30*4),30,4)    

test_that("create.life cannot get dimensions of dataset", {
	expect_that(create.life(lv.coefs = lv.coefs, family = "normal"), throws_error())
	})
	
test_that("create.life needs both X and X.coefs", {
	expect_that(create.life(lv.coefs = lv.coefs, X.coefs = X.coefs, family = "normal"), throws_error())
	})

sim.y <- create.life(true.lv, lv.coefs, X, X.coefs, family = "normal")
 
test_that("create.life default", {
	expect_that(dim(sim.y), equals(c(30,30)))
	})
 

library(mvabund)
data(antTraits)
traits <- as.matrix(cbind(1,antTraits$traits[,c(1,2,5)]))[1:30,]
traits.coefs <- cbind(matrix(rnorm((ncol(X)+1)*4),ncol=4),1)

test_that("create.life needs both traits and which.traits", {
	expect_that(create.life(true.lv, lv.coefs, X, X.coefs, traits = traits, family = "normal"), throws_error())
	expect_that(create.life(true.lv, lv.coefs, X, X.coefs, traits.coefs = traits.coefs, family = "normal"), throws_error())
	})

test_that("create.life will generate new coefs if traits supplied", {
	expect_that(create.life(true.lv, lv.coefs, X, X.coefs, traits.coefs = traits.coefs, traits = traits, family = "normal"), gives_warning())
	})
	
test_that("create.life has correct number of familes", {
	expect_that(create.life(true.lv, lv.coefs, X, X.coefs, family = c("poisson","normal")), throws_error())
	})

test_that("create.life needs trial size right", {
	expect_that(create.life(true.lv, lv.coefs, X, X.coefs, family = "binomial", trial.size = c(1,2)), throws_error())
	})
	
test_that("create.life needs cutoffs for ordinal", {
	expect_that(create.life(true.lv, lv.coefs, X, X.coefs, family = "ordinal", trial.size = c(1,2)), throws_error())
	})
	

test_that("create.life needs row.params correct", {
	expect_that(create.life(true.lv, lv.coefs, X, X.coefs, row.eff = "fixed", family = "poisson"), throws_error())
	expect_that(create.life(true.lv, lv.coefs, X, X.coefs, row.eff = "random", row.params = c(1,2), family = "poisson"), throws_error())
	})
	
	
		
################
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
X <- scale(spider$x)
n <- nrow(y); p <- ncol(y); 


fit1 <- boral(y, family = "negative.binomial", num.lv = 2, n.burnin = 10, n.iteration = 100, n.thin = 1, row.eff = "fixed", save.model = TRUE, calc.ics = TRUE)
getresiduals <- ds.residuals(fit1)

test_that("residuals default", {
	expect_that(length(getresiduals), equals(2))
	expect_that(dim(getresiduals$residuals), equals(c(n,p)))
	})

	
############
fit1 <- boral(y, family = "negative.binomial", num.lv = 2, n.burnin = 10, n.iteration = 100, n.thin = 1, row.eff = "fixed", save.model = TRUE, calc.ics = TRUE)
getfitted <- fitted(fit1)

test_that("fitted default", {
	expect_that(length(getfitted), equals(2))
	expect_that(dim(getfitted$out), equals(c(n,p)))
	})

	
#################	
fit1 <- boral(y, family = "negative.binomial", num.lv = 2, n.burnin = 10, n.iteration = 100, n.thin = 1, row.eff = "fixed", save.model = TRUE, calc.ics = TRUE)
fit.mcmc <- mcmc(fit1$jags.model$BUGSoutput$sims.matrix)                

test_that("get.measures gets wrong family", {
	expect_that(get.measures(y, family = c("poisson","negative.binomial"), num.lv = fit1$num.lv, fit.mcmc = fit.mcmc, row.eff = "fixed"), throws_error())
	})

test_that("get.measures will not do more measures if num.lv = 0", {
	expect_that(get.measures(y, family = "negative.binomial", num.lv = 0, fit.mcmc = fit.mcmc, row.eff = "fixed", more.measures = TRUE), throws_error())
	})
	
	
#################	
fit1 <- boral(y, family = "negative.binomial", num.lv = 2, n.burnin = 10, n.iteration = 100, n.thin = 1, row.eff = "fixed", save.model = TRUE, calc.ics = TRUE)
fit.mcmc <- mcmc(fit1$jags.model$BUGSoutput$sims.matrix)                

test_that("get.more.measures won't do anything if num.lv = 0", {
	expect_that(get.more.measures(y, family = "negative.binomial", num.lv = 0, fit.mcmc = fit.mcmc, row.eff = "fixed"), throws_error())
	})

test_that("get.more.measures gets wrong family", {
	expect_that(get.more.measures(y, family = c("poisson","negative.binomial"), num.lv = fit1$num.lv, fit.mcmc = fit.mcmc, row.eff = "fixed"), throws_error())
	})

test_that("get.more.measures gets wrong # of lvs", {
	skip_on_cran()
	expect_that(get.more.measures(y, family = "negative.binomial", num.lv = 1, fit.mcmc = fit.mcmc, row.eff = "fixed"), throws_error())
	})


test_that("get.more.measures default", {
	skip_on_cran()
	do.measures <- get.more.measures(y, family = "negative.binomial", num.lv = fit1$num.lv, fit.mcmc = fit.mcmc, row.eff = "fixed")                
	expect_that(length(do.measures), equals(6))
	})
	
	
fit1 <- boral(y, family = "negative.binomial", num.lv = 2, n.burnin = 10, n.iteration = 100, n.thin = 1, row.eff = "fixed", save.model = TRUE, calc.ics = TRUE)
fit.mcmc <- mcmc(fit1$jags.model$BUGSoutput$sims.matrix)


test_that("get.more.measures default", {
	skip_on_cran()
	expect_that(length(get.more.measures(y, family = "negative.binomial", num.lv = fit1$num.lv, fit.mcmc = fit.mcmc, row.eff = "fixed")), equals(6))
	})		

		
###################
fit1 <- boral(y, X = X, family = "negative.binomial", num.lv = 2, row.eff = "fixed", n.burnin = 10, n.iteration = 100, n.thin = 1, calc.ics = TRUE)
  	
test_that("get.enviro.cor needs MCMC samples", {
	expect_that(get.enviro.cor(fit1), throws_error())
	})

fit1 <- boral(y, family = "negative.binomial", num.lv = 2, row.eff = "fixed", n.burnin = 10, n.iteration = 100, n.thin = 1, calc.ics = TRUE)
  	
test_that("get.enviro.cor needs MCMC samples for X", {
	expect_that(get.enviro.cor(fit1), throws_error())
	})

fit1 <- boral(y, X = X, family = "negative.binomial", num.lv = 2, n.burnin = 10, n.iteration = 100, n.thin = 1, calc.ics = TRUE, save.model = TRUE)
	
ecors <- get.enviro.cor(fit1)

test_that("get.enviro.cor gets dimensions right", {
	expect_that(dim(ecors$cor), equals(c(p,p)))
	expect_that(dim(ecors$sig.cor), equals(c(p,p)))
	expect_that(dim(ecors$cov), equals(c(p,p)))
	})
	

#################
fit1 <- boral(y, X = X, family = "negative.binomial", num.lv = 2, row.eff = "fixed", n.burnin = 10, n.iteration = 100, n.thin = 1, calc.ics = TRUE)

test_that("get.residual.cor needs MCMC samples", {
	expect_that(get.residual.cor(fit1), throws_error())
	})

fit1 <- boral(y, X = X, family = "negative.binomial", row.eff = "fixed", n.burnin = 10, n.iteration = 100, n.thin = 1, calc.ics = TRUE, save.model = TRUE)

test_that("get.residual.cor needs MCMC samples of lv", {
	expect_that(get.residual.cor(fit1), throws_error())
	})

	
fit1 <- boral(y, X = X, family = "negative.binomial", num.lv = 2, row.eff = "fixed", n.burnin = 10, n.iteration = 100, n.thin = 1, calc.ics = TRUE, save.model = TRUE)
	
res.cors <- get.residual.cor(fit1)

test_that("get.residual.cor gets dimensions right", {
	expect_that(dim(res.cors$cor), equals(c(p,p)))
	expect_that(dim(res.cors$sig.cor), equals(c(p,p)))
	expect_that(dim(res.cors$cov), equals(c(p,p)))
	})

				
#################
test_that("makejagsboralnullmodel needs both X and traits", {
	expect_that(make.jagsboralnullmodel(family = rep(c("poisson","negative.binomial"),length=p), num.X = 0, num.traits = 5, row.eff = "fixed", n = n, p = p), throws_error())
	})
	

test_that("makejagsboralnullmodel which.traits wrong length", {
	expect_that(make.jagsboralnullmodel(family = rep(c("poisson","negative.binomial"),length=p), num.X = 2, num.traits = 5, which.traits = list(a = c(1,3)), row.eff = "fixed", n = n, p = p), throws_error())
	})

test_that("makejagsboralnullmodel has wrong trial size", {
	expect_that(make.jagsboralnullmodel(family = "binomial", num.X = 2, trial.size = c(1,2), row.eff = "fixed", n = n, p = p), throws_error())
	})

test_that("makejagsboralnullmodel has wrong family length", {
	expect_that(make.jagsboralnullmodel(family = c("poisson","binomial"), num.X = 2, trial.size = c(1,2), row.eff = "fixed", n = n, p = p), throws_error())
	})

	
#################
test_that("makejagsboralmodel needs both X and traits", {
	expect_that(make.jagsboralmodel(family = rep(c("poisson","negative.binomial"),length=p), num.X = 0, num.traits = 5, row.eff = "fixed", n = n, p = p), throws_error())
	})
	

test_that("makejagsboralmodel which.traits elements are wrong length", {
	expect_that(make.jagsboralmodel(family = rep(c("poisson","negative.binomial"),length=p), num.X = 2, num.traits = 1, which.traits = list(a = c(1,3)), row.eff = "fixed", n = n, p = p), throws_error())
	})

test_that("makejagsboralmodel has wrong trial size", {
	expect_that(make.jagsboralmodel(family = "binomial", num.X = 2, trial.size = c(1,2), row.eff = "fixed", n = n, p = p), throws_error())
	})

test_that("makejagsboralmodel has wrong family length", {
	expect_that(make.jagsboralmodel(family = c("poisson","binomial"), num.X = 2, trial.size = c(1,2), row.eff = "fixed", n = n, p = p), throws_error())
	})
	

################
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
X <- scale(spider$x)
n <- nrow(y); p <- ncol(y); 

fit1 <- boral(y, X = X, family = "negative.binomial", num.lv = 2, row.eff = "fixed", n.burnin = 10, n.iteration = 100, n.thin = 1, calc.ics = TRUE)

fit2 <- manyglm(mvabund(y) ~ X)

test_that("simulate function needs boral", {
	expect_that(simulate(fit2), throws_error())
	})

trysim <- simulate(fit1,nsim=10)

test_that("simulate function", {
	expect_that(dim(trysim), equals(c(nrow(y),ncol(y),10)))
	})

	