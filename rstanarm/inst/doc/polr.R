## ---- SETTINGS-knitr, include=FALSE--------------------------------------
stopifnot(require(knitr))
opts_chunk$set(
  comment=NA, message = FALSE, warning = FALSE,
  fig.align='center', fig.width = 7, fig.height = 3
)

## ---- SETTINGS-gg, include=FALSE-----------------------------------------
library(ggplot2)
thm_els <- theme(axis.text.y = element_blank(), 
                 legend.position = "none",
                 legend.background = element_rect(fill = "gray"),
                 legend.text = element_text(size = 7))
theme_set(theme_classic() %+replace% thm_els)

## ---- SETTINGS-rstan, include=FALSE--------------------------------------
ITER <- 500L
CHAINS <- 2L
CORES <- 1L
SEED <- 12345

## ---- SETTINGS-loo, include=FALSE----------------------------------------
loo.cores <- if (exists("CORES")) CORES else 1L
options(loo.cores = loo.cores)

## ----polr-tobgp-mass-----------------------------------------------------
library(MASS)
print(polr(tobgp ~ agegp + alcgp, data = esoph), digits = 1)

## ----polr-tobgp-mcmc, results="hide"-------------------------------------
library(rstanarm)
post0 <- stan_polr(tobgp ~ agegp + alcgp, data = esoph, 
                   prior = R2(0.25), prior_counts = dirichlet(1),
                   chains = CHAINS, cores = CORES, seed = SEED, iter = 200)

## ------------------------------------------------------------------------
print(post0, digits = 1)

## ---- polr-tobgp-cutpoints, echo=FALSE-----------------------------------
zeta_medians <- round(apply(rstan::extract(post0$stanfit, pars = "zeta")[[1]], 
                            2, median), digits = 2)

## ----polr-birthwt-recodes------------------------------------------------
data("birthwt", package = "MASS")
birthwt$race <- factor(birthwt$race, levels = 1:3, 
                       labels = c("white", "black", "other"))
birthwt$bwt <- birthwt$bwt / 1000 # convert from grams to kilograms
birthwt$low <- factor(birthwt$low, levels = 0:1, labels = c("no", "yes"))

## ----polr-stan_lm, results="hide"----------------------------------------
post1 <- stan_lm(-bwt ~ smoke + age + race + ptl + ht + ftv,
                 data = birthwt, prior = R2(0.5), 
                 chains = CHAINS, cores = CORES, seed = SEED)

## ------------------------------------------------------------------------
print(post1)

## ----polr-birthwt-mcmc, results="hide"-----------------------------------
post2 <- stan_polr(low ~ smoke + age + race + ptl + ht + ftv, data = birthwt,
                   prior = R2(0.5), prior_counts = dirichlet(c(1,1)), 
                   method = "probit", 
                   chains = CHAINS, cores = CORES, seed = SEED)

## ---- polr-loo-plot------------------------------------------------------
plot(loo(post2))

## ----polr-birthwt-comparison---------------------------------------------
round(cbind(Linear = coef(post1), Ordinal = coef(post2), 
            Rescaled = coef(post1) / sigma(post1)), 3)


