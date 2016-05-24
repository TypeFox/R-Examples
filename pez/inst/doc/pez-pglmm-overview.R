## ----include=FALSE--------------------
require(pez)
options(width=40)

## ----tidy=TRUE------------------------
nspp <- 15
nsite <- 10
env <- 1:nsite
env <- as.numeric(scale(env))

## ----tidy=TRUE------------------------
require(pez)
require(ape)
phy <- rcoal(n=nspp)
Vphy <- vcv(phy)
Vphy <- Vphy/(det(Vphy)^(1/nspp))

## ----tidy=TRUE------------------------
iD <- t(chol(Vphy))
intercept <- iD %*% rnorm(nspp)
slope <- iD %*% rnorm(nspp)

## ----tidy=TRUE------------------------
prob <- rep(intercept, each=nsite)
prob <- prob + rep(slope, each=nsite) * rep(env, nspp)
prob <- prob + rnorm(nspp*nsite)
pres <- rbinom(length(prob), size=1, prob=exp(prob)/(1+exp(prob)))

## ----tidy=TRUE------------------------
site <- factor(rep(1:nsite, each=nspp))
species <- factor(rep(1:nspp, nsite))
env <- rep(env, nspp)

## ----tidy=TRUE------------------------
r.intercept.spp.indep <- list(1, sp = species, covar = diag(nspp))
r.intercept.spp.phy <- list(1, sp = species, covar = Vphy)
r.slope.spp.indep <- list(env, sp = species, covar = diag(nspp))
r.slope.spp.phy <- list(env, sp = species, covar = Vphy)   
r.site <- list(1, site = site, covar = diag(nsite))
rnd.effects <- list(r.intercept.spp.indep, r.intercept.spp.phy, r.slope.spp.indep, r.slope.spp.phy, r.site)

## ----tidy=TRUE------------------------
model <- communityPGLMM(pres ~ env, family = "binomial", sp = species, site = site, random.effects = rnd.effects, REML = TRUE, verbose = FALSE)
communityPGLMM.binary.LRT(model, re.number = 1)
communityPGLMM.binary.LRT(model, re.number = 2)

