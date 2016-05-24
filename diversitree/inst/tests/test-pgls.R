source("helper-diversitree.R")

context("PGLS")

## Simulated tree and traits:
set.seed(1)
phy <- tree.bd(pars=c(1,0), max.taxa=100)
ty <- sim.character(phy, 1) + 3
ta <- sim.character(phy, 1)
tb <- sim.character(phy, 1)

data <- data.frame(a=ta, b=tb, y=ty, row.names=names(ty))

## Fit the model using standard approaches:
## 1. caper::pgls
cdata <- comparative.data(phy, cbind(data, species=phy$tip.label),
                          "species")
fit.caper.ya  <- pgls(y ~ a, cdata)
fit.caper.yab <- pgls(y ~ a + b, cdata)

## 2. nlme::gls
fit.gls.ya  <- gls(y ~ a,     data, corBrownian(1, phy), method="ML")
fit.gls.yab <- gls(y ~ a + b, data, corBrownian(1, phy), method="ML")

## This is a direct simplification of the code from Freckleton (2012);
## skipping missing data treatment, and returning the ML diffusion
## estimate.  It is this last part that we need here.
pgls.contrasts <- function(formula, data, phylo) {
  cdata <- data.frame(apply(data, 2, pic, phylo))
  fm <- lm(formula, cdata, y=TRUE)
  vars <- pic(data[[1]], phylo, var.contrasts=TRUE)[,2]
  V <- diversitree:::pgls.root.var.bm(phylo)
  u <- residuals(fm) 
  n <- length(u)
  sigma2 <- sum(u^2) / (n + 1)
  logLik <- -0.5 *( (n + 1) * log( 2 * pi * sigma2) + sum(log( vars ))
                   +  log(V) + sum(u^2) / sigma2)
  list(model=fm, logLik=logLik, sigma2=sigma2)
}

fit.contrasts.ya  <- pgls.contrasts(y ~ a     - 1, data, phy)
fit.contrasts.yab <- pgls.contrasts(y ~ a + b - 1, data, phy)

test_that("Reference implementations agree", {
  expect_that(fit.contrasts.ya$logLik,
              equals(as.numeric(logLik(fit.gls.ya))))
  expect_that(fit.contrasts.yab$logLik,
              equals(as.numeric(logLik(fit.gls.yab))))

  expect_that(coef(fit.contrasts.ya$model),
              equals(coef(fit.caper.ya)[-1]))
  expect_that(coef(fit.contrasts.yab$model),
              equals(coef(fit.caper.yab)[-1]))

  expect_that(coef(fit.gls.ya),  equals(coef(fit.caper.ya)))
  expect_that(coef(fit.gls.yab), equals(coef(fit.caper.yab)))
})

## From this, extract the full set of ML parameters (this requires
## interrogating both the caper and contrasts fits)
p.ya <- unname(c(coef(fit.caper.ya), fit.contrasts.ya$sigma2))
p.yab <- unname(c(coef(fit.caper.yab), fit.contrasts.yab$sigma2))
l.ya <- fit.contrasts.ya$logLik
l.yab <- fit.contrasts.yab$logLik

## Arbitrary offset:
p2.ya  <- p.ya  + runif(length(p.ya),  0, 0.1)
p2.yab <- p.yab + runif(length(p.yab), 0, 0.1)

## 4. VCV equation
lik.vcv.ya  <- make.pgls(phy, y ~ a,     data, control=list(method="vcv"))
lik.vcv.yab <- make.pgls(phy, y ~ a + b, data, control=list(method="vcv"))

test_that("VCV Likelihoods agree at ML points", {
  expect_that(lik.vcv.ya(p.ya),   equals(l.ya))
  expect_that(lik.vcv.yab(p.yab), equals(l.yab))
})

## 5. Contrasts
lik.con.ya  <- make.pgls(phy, y ~ a,     data, control=list(method="contrasts"))
lik.con.yab <- make.pgls(phy, y ~ a + b, data, control=list(method="contrasts"))

test_that("Contrasts Likelihoods agree at ML points", {
  expect_that(lik.con.ya(p.ya),   equals(l.ya))
  expect_that(lik.con.yab(p.yab), equals(l.yab))
})

test_that("Calculations agree at differing point", {
  expect_that(lik.con.ya(p2.ya),   equals(lik.vcv.ya(p2.ya)))
  expect_that(lik.con.yab(p2.yab), equals(lik.vcv.yab(p2.yab)))
})

## 6. Pruning
lik.pru.R.ya  <- make.pgls(phy, y ~ a,     data,
                           control=list(method="pruning", backend="R"))
lik.pru.R.yab <- make.pgls(phy, y ~ a + b, data,
                           control=list(method="pruning", backend="R"))
lik.pru.C.ya  <- make.pgls(phy, y ~ a,     data,
                           control=list(method="pruning", backend="C"))
lik.pru.C.yab <- make.pgls(phy, y ~ a + b, data,
                           control=list(method="pruning", backend="C"))

test_that("Pruning likelihoods agree at ML points", {
  expect_that(lik.pru.R.ya(p.ya),   equals(l.ya))
  expect_that(lik.pru.R.yab(p.yab), equals(l.yab))
  expect_that(lik.pru.C.ya(p.ya),   equals(l.ya))
  expect_that(lik.pru.C.yab(p.yab), equals(l.yab))
})

test_that("Calculations agree at differing point", {
  expect_that(lik.pru.R.ya(p2.ya),   equals(lik.vcv.ya(p2.ya)))
  expect_that(lik.pru.R.yab(p2.yab), equals(lik.vcv.yab(p2.yab)))
  expect_that(lik.pru.C.ya(p2.ya),   equals(lik.vcv.ya(p2.ya)))
  expect_that(lik.pru.C.yab(p2.yab), equals(lik.vcv.yab(p2.yab)))
})

test_that("Fitted values are correct", {
  expect_that(fitted(lik.vcv.ya, p.ya),     equals(fitted(fit.caper.ya)))
  expect_that(fitted(lik.con.ya, p.ya),     equals(fitted(fit.caper.ya)))
  expect_that(fitted(lik.pru.R.ya, p.ya),   equals(fitted(fit.caper.ya)))
  expect_that(fitted(lik.pru.C.ya, p.ya),   equals(fitted(fit.caper.ya)))
  expect_that(fitted(lik.vcv.yab, p.yab),   equals(fitted(fit.caper.yab)))
  expect_that(fitted(lik.con.yab, p.yab),   equals(fitted(fit.caper.yab)))
  expect_that(fitted(lik.pru.R.yab, p.yab), equals(fitted(fit.caper.yab)))
  expect_that(fitted(lik.pru.C.yab, p.yab), equals(fitted(fit.caper.yab)))
})

test_that("Residuals values are correct", {
  expect_that(resid(lik.vcv.ya, p.ya),     equals(resid(fit.caper.ya)))
  expect_that(resid(lik.con.ya, p.ya),     equals(resid(fit.caper.ya)))
  expect_that(resid(lik.pru.R.ya, p.ya),   equals(resid(fit.caper.ya)))
  expect_that(resid(lik.pru.C.ya, p.ya),   equals(resid(fit.caper.ya)))
  expect_that(resid(lik.vcv.yab, p.yab),   equals(resid(fit.caper.yab)))
  expect_that(resid(lik.con.yab, p.yab),   equals(resid(fit.caper.yab)))
  expect_that(resid(lik.pru.R.yab, p.yab), equals(resid(fit.caper.yab)))
  expect_that(resid(lik.pru.C.yab, p.yab), equals(resid(fit.caper.yab)))
})

fit.ya  <- find.mle(lik.con.ya, p.ya)
fit.yab <- find.mle(lik.con.yab, p.yab)

test_that("Fitted values are correct (ML search object)", {
  expect_that(fitted(fit.ya),
              is_identical_to(fitted(lik.con.ya, coef(fit.ya))))
  expect_that(fitted(fit.yab),
              is_identical_to(fitted(lik.con.yab, coef(fit.yab))))
})

test_that("Residuals values are correct (ML search object)", {
  expect_that(resid(fit.ya),
              is_identical_to(resid(lik.con.ya, coef(fit.ya))))
  expect_that(resid(fit.yab),
              is_identical_to(resid(lik.con.yab, coef(fit.yab))))
})

set.seed(1)
samples.ya  <- mcmc(lik.con.ya, coef(fit.ya), 100, w=1, print.every=0)
samples.yab <- mcmc(lik.con.yab, coef(fit.yab), 100, w=1, print.every=0)

test_that("Fitted values are correct (MCMC object)", {
  m2l <- diversitree:::matrix.to.list
  f <- function(p, lik)
    do.call(cbind, lapply(m2l(p), function(pi) fitted(lik, pi)))
  
  expect_that(fitted(samples.ya),
              is_identical_to(f(coef(samples.ya), lik.con.ya)))
  expect_that(fitted(samples.yab),
              is_identical_to(f(coef(samples.yab), lik.con.yab)))

  ## Argument passing:
  burnin <- 10
  thin <- 2
  sample <- 10
  set.seed(1)
  cmp <- f(coef(samples.ya, burnin=burnin, thin=thin, sample=sample),
           lik.con.ya)
  set.seed(1)
  expect_that(fitted(samples.ya, burnin=burnin, thin=thin, sample=sample),
              is_identical_to(cmp))
})

test_that("Residuals are correct (MCMC object)", {
  m2l <- diversitree:::matrix.to.list
  f <- function(p, lik)
    do.call(cbind, lapply(m2l(p), function(pi) resid(lik, pi)))
  
  expect_that(resid(samples.ya),
              is_identical_to(f(coef(samples.ya), lik.con.ya)))
  expect_that(resid(samples.yab),
              is_identical_to(f(coef(samples.yab), lik.con.yab)))

  ## Argument passing:
  burnin <- 10
  thin <- 2
  sample <- 10
  set.seed(1)
  cmp <- f(coef(samples.ya, burnin=burnin, thin=thin, sample=sample),
           lik.con.ya)
  set.seed(1)
  expect_that(resid(samples.ya, burnin=burnin, thin=thin, sample=sample),
              is_identical_to(cmp))
})
