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

## ---- count-roaches-mcmc, results="hide"---------------------------------
library(rstanarm)
data(roaches)

# Rescale
roaches$roach1 <- roaches$roach1 / 100
# Estimate original model
glm1 <- glm(y ~ roach1 + treatment + senior, offset = log(exposure2), 
            data = roaches, family = poisson)
# Estimate Bayesian version with stan_glm
stan_glm1 <- stan_glm(y ~ roach1 + treatment + senior, offset = log(exposure2),
                      data = roaches, family = poisson, 
                      prior = normal(0,2.5), prior_intercept = normal(0,5),
                      chains = CHAINS, cores = CORES, seed = SEED)

## ---- count-roaches-comparison-------------------------------------------
round(rbind(glm = coef(glm1), stan_glm = coef(stan_glm1)), digits = 2)
round(rbind(glm = summary(glm1)$coefficients[, "Std. Error"], 
            stan_glm = se(stan_glm1)), digits = 3)

## ---- count-roaches-posterior_predict------------------------------------
yrep <- posterior_predict(stan_glm1)

## ---- count-roaches-plot-pp_check1, fig.height=3, fig.width=4------------
prop_zero <- function(y) mean(y == 0)
(prop_zero_test1 <- pp_check(stan_glm1, check = "test", test = "prop_zero"))

## ---- count-roaches-negbin, results="hide"-------------------------------
stan_glm2 <- update(stan_glm1, family = neg_binomial_2) 

## ---- count-roaches-plot-pp_check2, fig.height=3, fig.width=8------------
library(gridExtra)
prop_zero_test2 <- pp_check(stan_glm2, check = "test", test = "prop_zero")
# Show graphs for Poisson and negative binomial side by side
grid.arrange(prop_zero_test1 + ggtitle("Poisson"), 
             prop_zero_test2 + ggtitle("Negative Binomial"), 
             ncol = 2)

## ---- count-roaches-loo--------------------------------------------------
loo1 <- loo(stan_glm1)
loo2 <- loo(stan_glm2)
compare(loo1, loo2)

