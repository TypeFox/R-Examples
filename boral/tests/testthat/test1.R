context("test boral.jags.R")

################
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
X <- scale(spider$x)
n <- nrow(y); p <- ncol(y); 

fit1 <- boral(y, X = X, family = "negative.binomial", num.lv = 2, row.eff = "fixed", n.burnin = 10, n.iteration = 100, n.thin = 1, calc.ics = TRUE)

test_that("boral.default attributes", {
	expect_that(fit1, is_a("boral"))
	expect_that(length(fit1$family), equals(p))
	expect_that(dim(fit1$y), equals(c(n,p)))
	expect_that(fit1$num.X, equals(ncol(X)))
	expect_that(fit1$num.lv, equals(2))
	expect_that(fit1$hypparams, equals(c(100, 20, 100, 50)))
	expect_that(c(fit1$n.burnin, fit1$n.thin, fit1$n.iteration), equals(c(10, 1, 100)))
	expect_that(length(fit1$ics), equals(6))	
	expect_that(fit1$row.eff, matches("fixed"))	
	expect_that(fit1$ssvs.index, equals(rep(-1,ncol(X))))
	expect_that(fit1$num.ord.levels, equals(0))
	})


test_that("too few or many hypparams", {
	expect_that(boral(y, X = X, family = "negative.binomial", num.lv = 2, row.eff = "fixed", n.burnin = 10, n.iteration = 100, n.thin = 1, hypparams = c(100,100)), throws_error())
	})

test_that("too few or many lvs", {
	expect_that(boral(y, X = X, family = "negative.binomial", num.lv = 10, row.eff = "fixed", n.burnin = 10, n.iteration = 100, n.thin = 1), gives_warning())
	expect_that(boral(y, X = X, family = "negative.binomial", num.lv = 1, row.eff = "fixed", n.burnin = 10, n.iteration = 100, n.thin = 1), gives_warning())
	})

test_that("too length of families not correct", {
	expect_that(boral(y, X = X, family = c("poisson","negative.binomial"), row.eff = "fixed", n.burnin = 10, n.iteration = 100, n.thin = 1), throws_error())
	})

test_that("row effects wrong", {
	expect_that(boral(y, X = X, family = "poisson", row.eff = "blah", n.burnin = 10, n.iteration = 100, n.thin = 1), throws_error())
	})

test_that("ssvs.index incorrect length or values wrong", {
	expect_that(boral(y, X = X, family = "negative.binomial", ssvs.index = c(1,2), n.burnin = 10, n.iteration = 100, n.thin = 1), throws_error())
	expect_that(boral(y, X = X, family = "negative.binomial", ssvs.index = c(-2), n.burnin = 10, n.iteration = 100, n.thin = 1), throws_error())
	})
	
	
test_that("trial.size incorrect", {
	expect_that(boral(y, X = X, family = "binomial", trial.size = c(1,2,3), n.burnin = 10, n.iteration = 100, n.thin = 1), throws_error())
	})

test_that("do.fit = FALSE default", {
	expect_that(boral(y, X = X, family = "negative.binomial", num.lv = 2, row.eff = "fixed", n.burnin = 10, n.iteration = 100, n.thin = 1, do.fit = FALSE), prints_text("JAGS model file created only. Thank you, come again!"))
	})

	
###############
test_that("lvsplot default", {
	expect_that(lvsplot(fit1,ind.spp=5), prints_text("Only the first 5 `most important' latent variable coefficients included in biplot"))
	})

	
fit1 <- boral(y, X = X, family = "negative.binomial", num.lv = 0, row.eff = "fixed", n.burnin = 10, n.iteration = 100, n.thin = 1)
test_that("lvsplot no lvs", {
	expect_that(lvsplot(fit1), throws_error())
	})
	
fit1 <- boral(y, family = "negative.binomial", num.lv = 3, row.eff = "fixed", n.burnin = 10, n.iteration = 100, n.thin = 1, calc.ics = FALSE)
test_that("lvsplot cannot exceed 2 lvs", {
	expect_that(lvsplot(fit1), throws_error())
	})	
	
###############
fit1 <- boral(y, X = X, family = "negative.binomial", num.lv = 2, row.eff = "fixed", n.burnin = 10, n.iteration = 100, n.thin = 1, calc.ics = TRUE)
sumfit1 <- summary(fit1)	

test_that("summary default", {
	expect_that(sumfit1, is_a("summary.boral"))
	expect_that(sumfit1$num.ord.levels, equals(0))	
	expect_that(sumfit1$est, matches("median"))	
	expect_that(length(sumfit1$ics), equals(6))	
	expect_that(sumfit1$ssvs.index, equals(rep(-1,ncol(X))))
	expect_that(sumfit1$num.ord.levels, equals(0))
	})

	
###############
X <- cbind(1,scale(spider$x))
n <- nrow(y); p <- ncol(y); 

test_that("boral does not permit intercept", {
	expect_that(boral(y, X = X, family = "negative.binomial", num.lv = 2, row.eff = "fixed", n.burnin = 10, n.iteration = 100, n.thin = 1, calc.ics = FALSE), throws_error())
	})
	
	
###############
data(antTraits)
y <- antTraits$abun
X <- as.matrix(scale(antTraits$env))
traits <- as.matrix(cbind(1,antTraits$traits[,c(1,2,5)]))
which.traits <- vector("list",ncol(X)+1)
for(i in 1:length(which.traits)) which.traits[[i]] <- 1:ncol(traits)
    
test_that("boral needs X supplied if traits supplied", {
	expect_that(boral(y, traits = traits, which.traits = which.traits, family = "negative.binomial", num.lv = 2, n.burnin = 10, n.iteration = 100, n.thin = 1, calc.ics = FALSE), throws_error())
	})
	
test_that("supply which.traits if traits supplied", {
	expect_that(boral(y, X = X, traits = traits, which.traits = NULL, family = "negative.binomial", num.lv = 2, n.burnin = 10, n.iteration = 100, n.thin = 1, calc.ics = FALSE), throws_error())
	})
	

which.traits <- vector("list",2)
for(i in 1:length(which.traits)) which.traits[[i]] <- 1:ncol(traits)

test_that("which.traits is of wrong life", {
	expect_that(boral(y, X = X, traits = traits, which.traits = NULL, family = "negative.binomial", num.lv = 2, n.burnin = 10, n.iteration = 100, n.thin = 1, calc.ics = FALSE), throws_error())
	})
	
###############
y <- matrix(sample(1:5,50,replace=T),10)

fit1 <- boral(y, family = "ordinal", num.lv = 0, row.eff = "fixed", n.burnin = 10, n.iteration = 100, n.thin = 1, calc.ics = FALSE)

test_that("cannot do plot.boral for ordinal", {
	expect_that(plot(fit1), throws_error())
	})

