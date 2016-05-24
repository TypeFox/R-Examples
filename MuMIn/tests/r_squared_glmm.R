if(length(find.package(c("lme4", "nlme"), quiet = TRUE)) == 2) {

library(lme4)
library(nlme)
library(MuMIn) ## needs  version >= 1.9.22

## Paul Johnson's code from MEE ms appendix:
R2PJ <-
function(model) { 
  X <- model.matrix(model)
  n <- nrow(X)
  Beta <- fixef(model)
  Sf <- var(X %*% Beta)
  Sigma.list <- VarCorr(model)
  Sl <- sum(
      sapply(Sigma.list, function(Sigma) {
          Xstar <- X[,rownames(Sigma)]
          sum(diag(Xstar %*% Sigma %*% t(Xstar))) / n
        }))
  Se <- attr(Sigma.list, "sc")^2
  Sd <- 0
  total.var <- Sf + Sl + Se + Sd
  return(c(R2m = Sf / total.var, R2c = (Sf + Sl) / total.var))
}

R2PJpois <-
function(model, ovdVarName) { 
  X <- model.matrix(model)
  n <- nrow(X)
  Beta <- fixef(model)
  Sf <- var(X %*% Beta)
  Sigma.list <- VarCorr(model)
  Sl <- sum(sapply(Sigma.list[names(Sigma.list) != ovdVarName], function(Sigma) {
          Xstar <-X[,rownames(Sigma)]
          sum(diag(Xstar %*% Sigma %*% t(Xstar)))/n
        }))
  Se <- Sigma.list[[ovdVarName]][1L]
  Sd <- log(1 + 1 / exp(mean(X %*% Beta))) ### See footnote ###
  (total.var <- c(Sf + Sl + Se + Sd))
  return (c(R2m = c(Sf / total.var), R2c = c((Sf + Sl) / total.var)))
}


### tests:
orangemod.ri <- lmer(circumference ~ age +(1 | Tree), data = Orange)
orangemod.rs <- lmer(circumference ~ age +(age | Tree), data = Orange)

stopifnot(
all.equal(r.squaredGLMM(orangemod.ri), R2PJ(orangemod.ri))
)

#r.squaredGLMM(orangemod.rs)
#R2PJ(orangemod.rs)

stopifnot(
all.equal(r.squaredGLMM(orangemod.rs), R2PJ(orangemod.rs))
)

set.seed(1)
Orange$rand <- rnorm(nrow(Orange))
subs <- 2:31
model1 <- lmer(circumference ~ scale(age) + (I(rand + 1) | Tree), data = Orange, subset = subs, REML = TRUE)
model1a <- lme(circumference ~ scale(age), ~ I(rand + 1) | Tree, data = Orange, subset = subs, method = "REML")
model1b <- lme(circumference ~ scale(age), list(~ I(rand + 1) | Tree), data = Orange, subset = subs, method = "REML")
rm(subs)

stopifnot(
all.equal(r.squaredGLMM(model1), r.squaredGLMM(model1a), tolerance = .001)
)

#fixef(model1)
#fixef(model1a)

model.frame(model1a, TRUE)
model.matrix(model1a, TRUE)

###########################

if(require(glmmML)) {

### tests models evaluated in a strange environment
#with(env <- new.env(), {
#
#	n <- 200
#	id <- factor(rep(1:20, rep(5, 20)))
#	y <- rbinom(n, prob = rep(runif(20), rep(5, 20)), size = 1)
#	x <- rnorm(n)
#	dat <- data.frame(y = y, x = x, id = id)
#	
#	subs <- 1:110
#	model2 <- glmmML(y ~ x, data = dat, cluster = id, subset = subs, x = TRUE)
#	model2a <- glmer(y ~ x + (1 | id), data = dat, subset = subs, family = family(model2))
#
#	model2b <- glmer(y ~ x + (1 | id), data = dat, subset = subs, family = binomial("logit"))
#	
#	
#	dat$..obslev <- gl(200, 1)
#	
#	x <- glmer(formula = y ~ x + (1 | id) + (1 | ..obslev), data = dat, family = binomial("logit"))
#	
#	model2b <- glmer(y ~ x + (1 | id) + (1 | gl(n, 1)), data = dat,  family = binomial("logit"),
#					 subset = subs)
#	
#	
#	model2b <- glmer(y ~ x + (1 | id) + (gl(200, 1), data = dat,  family = binomial("logit"))
#
#	x <- model2b
#
#})
#
#r.squaredGLMM(env$model2a)
#r.squaredGLMM(model2)

#stopifnot(
#all.equal(r.squaredGLMM(env$model2), r.squaredGLMM(env$model2a), tolerance = 0.005)
#)

} ## if(require(glmmML))
### Poisson:

# the value DUMMY=0 and the high counts DUMMY=1:
data(grouseticks, package = "lme4")
grouseticks <- merge(grouseticks, aggregate(
	list(medianTicks = grouseticks$TICKS), grouseticks['LOCATION'], median))
grouseticks$DUMMY <- as.integer(grouseticks$TICKS > grouseticks$medianTicks)
grouseticks$medianTicks <- NULL

tickmod.ri  <- glmer(TICKS ~ HEIGHT + DUMMY+YEAR + (1|INDEX) + (1|BROOD) + (1|LOCATION),
	family = poisson, data = grouseticks)
	
tickmod.ri0  <- glmer(update(formula(tickmod.ri), .~. -(1|INDEX)),
	family = poisson, data = grouseticks)

#x1  <- glm(TICKS ~ HEIGHT+DUMMY+YEAR, family = quasipoisson, data = grouseticks)
#
#x2 <- glmer(TICKS~HEIGHT+DUMMY+YEAR + (1|INDEX), family = poisson, data = grouseticks)
#
#sqrt(VarCorr(x2)[[1]][1])
#
#(11.47402)^2
#
#?QAIC
#qpf <- quasipoisson()
#qpf$aic <- poisson()$aic
#qpf$family <- "qp"
#
#tickmod.ri0q  <- glmer(TICKS~HEIGHT+DUMMY+YEAR + (1|BROOD) + (1|LOCATION),
#	family = qpf, data = grouseticks)
#tickmod.ri0q


print(R2PJpois(tickmod.ri, 'INDEX'))
print(r.squaredGLMM(tickmod.ri))
#print(r.squaredGLMM(tickmod.ri0))

stopifnot(all.equal(R2PJpois(tickmod.ri, 'INDEX'), r.squaredGLMM(tickmod.ri)))

## results are slightly different depending on the ordering of INDEX
## try with (1 | sample(INDEX))
stopifnot(
		  all.equal(r.squaredGLMM(tickmod.ri), r.squaredGLMM(tickmod.ri0), tolerance = 0.001)
		  )

tickmod.rs  <- glmer(TICKS ~ HEIGHT + DUMMY + YEAR + (1|INDEX) + (1|BROOD) + (DUMMY|LOCATION),
      family = poisson, data = grouseticks)
tickmod.rs0  <-  glmer(TICKS ~ HEIGHT + DUMMY + YEAR + (1|BROOD) + (DUMMY|LOCATION),
						  family = poisson, data = grouseticks)

stopifnot(all.equal(R2PJpois(tickmod.rs, 'INDEX'), r.squaredGLMM(tickmod.rs)))


#tickmod3 <- glmmML(TICKS ~ HEIGHT + DUMMY + YEAR, cluster = INDEX, family = poisson, data = grouseticks, x = TRUE)
#tickmod3a  <- glmer(TICKS ~ HEIGHT + DUMMY + YEAR + (1|INDEX), family = poisson, data = grouseticks)
#tickmod4 <- glmmML(TICKS ~ HEIGHT + DUMMY + YEAR, cluster = LOCATION, family = poisson, data = grouseticks, x = TRUE)
#r.squaredGLMM(tickmod3)
#r.squaredGLMM(tickmod3a)
#r.squaredGLMM(tickmod.rs0)
}