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

## ----lm-clouds-ols-------------------------------------------------------
data("clouds", package = "HSAUR3")
ols <- lm(rainfall ~ seeding * (sne + cloudcover + prewetness + echomotion) +
            time, data = clouds)
round(coef(ols), 3)

## ----lm-clouds-mcmc, results='hide'--------------------------------------
library(rstanarm)
post <- stan_lm(rainfall ~ seeding * (sne + cloudcover + prewetness + 
                                        echomotion) + time, data = clouds,
                prior = R2(location = 0.2), 
                chains = CHAINS, cores = CORES, seed = SEED)
post

## ---- echo=FALSE---------------------------------------------------------
print(post)

## ----lm-clouds-ate-plot, fig.height=3------------------------------------
clouds_cf <- clouds
clouds_cf$seeding[] <- "yes"
y1_rep <- posterior_predict(post, newdata = clouds_cf)
clouds_cf$seeding[] <- "no"
y0_rep <- posterior_predict(post, newdata = clouds_cf)
qplot(x = c(y1_rep - y0_rep), geom = "histogram", 
      ylab = NULL, xlab = "Estimated ATE")

## ----lm-clouds-simple, results="hide"------------------------------------
simple <- stan_glm(rainfall ~ seeding * (sne + cloudcover + prewetness + 
                                        echomotion) + time,
                   data = clouds, family = gaussian(), 
                   prior = cauchy(), prior_intercept = cauchy(),
                   chains = CHAINS, cores = CORES, seed = SEED)

## ----lm-clouds-loo, warning=TRUE-----------------------------------------
(loo_post <- loo(post))
(loo(simple))

## ----lm-clouds-plot-loo--------------------------------------------------
plot(loo_post, label_points = TRUE)

