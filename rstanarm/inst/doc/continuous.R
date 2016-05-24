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

## ---- continuous-kidiq-mcmc,results="hide"-------------------------------
library(rstanarm)
data(kidiq)
post1 <- stan_glm(kid_score ~ mom_hs, data = kidiq, family = gaussian(link = "identity"), 
                  chains = CHAINS, cores = CORES, seed = SEED, iter = ITER)
post2 <- update(post1, formula = . ~ mom_iq)
post3 <- update(post1, formula = . ~ mom_hs + mom_iq)
(post4 <- update(post1, formula = . ~ mom_hs * mom_iq))

## ---- continuous-kidiq-print, echo=FALSE---------------------------------
print(post4)

## ---- continuous-kidiq-plot1a, fig.width = 5-----------------------------
base <- ggplot(kidiq, aes(x = mom_hs, y = kid_score)) + 
  geom_point(size = 1, position = position_jitter(height = 0.05, width = 0.1)) + 
  scale_x_continuous(breaks = c(0,1), labels = c("No HS", "HS"))
  
base + geom_abline(intercept = coef(post1)[1], slope = coef(post1)[2], 
                   color = "skyblue4", size = 1)

## ---- continuous-kidiq-plot1b, fig.width = 5-----------------------------
draws <- as.data.frame(post1)
colnames(draws)[1:2] <- c("a", "b")

base + 
  geom_abline(data = draws, aes(intercept = a, slope = b), 
              color = "skyblue", size = 0.2, alpha = 0.25) + 
  geom_abline(intercept = coef(post1)[1], slope = coef(post1)[2], 
              color = "skyblue4", size = 1)

## ---- continuous-kidiq-plot2, fig.width = 5------------------------------
draws <- as.data.frame(as.matrix(post2))
colnames(draws)[1:2] <- c("a", "b")
ggplot(kidiq, aes(x = mom_iq, y = kid_score)) + 
  geom_point(size = 1) +
  geom_abline(data = draws, aes(intercept = a, slope = b), 
              color = "skyblue", size = 0.2, alpha = 0.25) + 
  geom_abline(intercept = coef(post2)[1], slope = coef(post2)[2], 
              color = "skyblue4", size = 1)

## ---- continuous-kidiq-plot3, fig.width = 5------------------------------
reg0 <- function(x, ests) cbind(1, 0, x) %*% ests 
reg1 <- function(x, ests) cbind(1, 1, x) %*% ests

args <- list(ests = coef(post3))
kidiq$clr <- factor(kidiq$mom_hs, labels = c("No HS", "HS"))
lgnd <- guide_legend(title = NULL)
base2 <- ggplot(kidiq, aes(x = mom_iq, fill = relevel(clr, ref = "HS"))) + 
  geom_point(aes(y = kid_score), shape = 21, stroke = .2, size = 1) + 
  guides(color = lgnd, fill = lgnd) + 
  theme(legend.position = "right")
base2 + 
  stat_function(fun = reg0, args = args, aes(color = "No HS"), size = 1.5) +
  stat_function(fun = reg1, args = args, aes(color = "HS"), size = 1.5)

## ---- continuous-kidiq-plot4, fig.width = 5------------------------------
reg0 <- function(x, ests) cbind(1, 0, x, 0 * x) %*% ests 
reg1 <- function(x, ests) cbind(1, 1, x, 1 * x) %*% ests
args <- list(ests = coef(post4))
base2 +
  stat_function(fun = reg0, args = args, aes(color = "No HS"), size = 1.5) + 
  stat_function(fun = reg1, args = args, aes(color = "HS"), size = 1.5)

## ---- continuous-kidiq-loo-----------------------------------------------
# Compare them with loo
loo1 <- loo(post1)
loo2 <- loo(post2)
loo3 <- loo(post3)
loo4 <- loo(post4)
compare(loo1, loo2, loo3, loo4)

## ---- continuous-kidiq-loo-2---------------------------------------------
compare(loo1, loo4)

## ---- continuous-kidiq-loo-3---------------------------------------------
compare(loo3, loo4)
compare(loo2, loo4)

## ---- continuous-kidiq-pp_check1, fig.width=5----------------------------
pp_check(post4, check = "distributions", overlay = FALSE, nreps = 5)

## ---- continuous-kidiq-pp_check2, fig.height = 3-------------------------
pp_check(post4, check = "test", test = "mean")

## ---- continuous-kidiq-pp_check3, fig.height=3---------------------------
pp_check(post4, check = "test", test = c("mean", "sd"))

## ---- continuous-kidiq-posterior_predict---------------------------------
IQ_SEQ <- seq(from = 75, to = 135, by = 5)
y_nohs <- posterior_predict(post4, newdata = data.frame(mom_hs = 0, mom_iq = IQ_SEQ))
y_hs <- posterior_predict(post4, newdata = data.frame(mom_hs = 1, mom_iq = IQ_SEQ))

## ---- continuous-kidiq-plot-predict, fig.width=9-------------------------
par(mfrow = c(1:2), mar = c(5,4,2,1))
boxplot(y_hs, axes = FALSE, outline = FALSE, ylim = c(10,170),
        xlab = "Mom IQ", ylab = "Predicted Kid IQ", main = "Mom HS")
axis(1, at = 1:ncol(y_hs), labels = IQ_SEQ, las = 3)
axis(2, las = 1)
boxplot(y_nohs, outline = FALSE, col = "red", axes = FALSE, ylim = c(10,170),
        xlab = "Mom IQ", ylab = NULL, main = "Mom No HS")
axis(1, at = 1:ncol(y_hs), labels = IQ_SEQ, las = 3)

## ---- continuous-kidiq-validation, eval=FALSE, include=FALSE-------------
#  # # External Validation
#  # source(paste0(ROOT, "ARM/Ch.3/kids_before1987.data.R"),
#  #        local = kidiq, verbose = FALSE)
#  # source(paste0(ROOT, "ARM/Ch.3/kids_after1987.data.R"),
#  #        local = kidiq, verbose = FALSE)
#  # post5 <- stan_lm(ppvt ~ hs + afqt, data = kidiq,
#  #                  prior = R2(location = 0.25, what = "mean"), seed = SEED)
#  # y_ev <- posterior_predict(post5, newdata =
#  #                           data.frame(hs = kidiq$hs_ev, afqt = kidiq$afqt_ev))
#  # par(mfrow = c(1,1))
#  # hist(-sweep(y_ev, 2, STATS = kidiq$ppvt_ev, FUN = "-"), prob = TRUE,
#  #      xlab = "Predictive Errors in ppvt", main = "", las = 2)

## ---- continuous-clotting-mle, results='hide'----------------------------
clotting <- data.frame(
    u = c(5,10,15,20,30,40,60,80,100),
    lot1 = c(118,58,42,35,27,25,21,19,18),
    lot2 = c(69,35,26,21,18,16,13,12,12))
summary(glm(lot1 ~ log(u), data = clotting, family = Gamma))
summary(glm(lot2 ~ log(u), data = clotting, family = Gamma))

## ---- continuous-clotting-mcmc, results="hide"---------------------------
clotting2 <- with(clotting, data.frame(
  log_plasma = rep(log(u), 2),
  clot_time = c(lot1, lot2),
  lot_id = factor(rep(c(1,2), each = length(u)))
))

fit <- stan_glm(clot_time ~ log_plasma * lot_id, data = clotting2, family = Gamma, 
                prior_intercept = normal(0,1), prior = normal(0,1),
                chains = CHAINS, cores = CORES, seed = SEED)

## ------------------------------------------------------------------------
print(fit, digits = 3)

