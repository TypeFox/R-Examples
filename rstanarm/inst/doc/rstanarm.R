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

## ----rstanarm-mle--------------------------------------------------------
data("womensrole", package = "HSAUR3")
womensrole$total <- womensrole$agree + womensrole$disagree
womensrole_glm_1 <- glm(cbind(agree, disagree) ~ education + gender,
                        data = womensrole, family = binomial(link = "logit"))
round(coef(summary(womensrole_glm_1)), 3)

## ----rstanarm-mcmc, results="hide"---------------------------------------
library(rstanarm)
womensrole_bglm_1 <- stan_glm(cbind(agree, disagree) ~ education + gender,
                              data = womensrole,
                              family = binomial(link = "logit"), 
                              prior = student_t(df = 7), 
                              prior_intercept = student_t(df = 7),
                              chains = CHAINS, cores = CORES, seed = SEED)
womensrole_bglm_1

## ---- echo=FALSE---------------------------------------------------------
print(womensrole_bglm_1)

## ----rstanarm-ci---------------------------------------------------------
ci95 <- posterior_interval(womensrole_bglm_1, prob = 0.95, pars = "education")
round(ci95, 2)

## ----rstanarm-methods----------------------------------------------------
cbind(Median = coef(womensrole_bglm_1), MAD_SD = se(womensrole_bglm_1))
summary(residuals(womensrole_bglm_1)) # not deviance residuals
cov2cor(vcov(womensrole_bglm_1))

## ----rstanarm-shinystan,eval=FALSE---------------------------------------
#  launch_shinystan(womensrole_bglm_1)

## ----rstanarm-posterior_predict------------------------------------------
y_rep <- posterior_predict(womensrole_bglm_1)
dim(y_rep)

## ----rstanarm-criticism-plot, fig.width=10, fig.cap="Posterior predictive boxplots vs. observed datapoints"----
par(mfrow = 1:2, mar = c(5,3.7,1,0) + 0.1, las = 3)
boxplot(sweep(y_rep[,womensrole$gender == "Male"], 2, STATS = 
               womensrole$total[womensrole$gender == "Male"], FUN = "/"), 
        axes = FALSE, main = "Male", pch = NA,
        xlab = "Years of Education", ylab = "Proportion of Agrees")
with(womensrole, axis(1, at = education[gender == "Male"] + 1, 
                      labels = 0:20))
axis(2, las = 1)
with(womensrole[womensrole$gender == "Male",], 
     points(education + 1,  agree / (agree + disagree), 
            pch = 16, col = "red"))
boxplot(sweep(y_rep[,womensrole$gender == "Female"], 2, STATS = 
          womensrole$total[womensrole$gender == "Female"], FUN = "/"), 
          axes = FALSE, main = "Female", pch = NA,
        xlab = "Years of Education", ylab = "")
with(womensrole, axis(1, at = education[gender == "Female"] + 1,
     labels = 0:20))
with(womensrole[womensrole$gender == "Female",], 
     points(education + 1,  agree / (agree + disagree), 
            pch = 16, col = "red"))

## ---- rstanarm-update, results="hide"------------------------------------
(womensrole_bglm_2 <- update(womensrole_bglm_1, formula. = . ~ . + I(education^2)))

## ---- echo=FALSE---------------------------------------------------------
print(womensrole_bglm_2)

## ----rstanarm-loo--------------------------------------------------------
loo_bglm_1 <- loo(womensrole_bglm_1)
loo_bglm_2 <- loo(womensrole_bglm_2)

## ----rstanarm-loo-plot, fig.width=7, fig.height=4------------------------
par(mfrow = 1:2, mar = c(5,3.8,1,0) + 0.1, las = 3)
plot(loo_bglm_1, label_points = TRUE)
plot(loo_bglm_2, label_points = TRUE)

## ---- rstanarm-loo-compare-----------------------------------------------
compare(loo_bglm_1, loo_bglm_2)

## ---- rstanarm-loo-print-------------------------------------------------
loo_bglm_1

## ---- rstanarm-posterior_predict-manipulate------------------------------
# note: in newdata we want agree and disgree to sum to the number of people we
# want to predict for. the values of agree and disagree don't matter so long as
# their sum is the desired number of trials. we need to explicitly imply the
# number of trials like this because our original data are aggregate. if we had
# bernoulli data then it would be a given we wanted to predict for single
# individuals.
newdata <- data.frame(agree = c(0,0), disagree = c(100,100), education = c(12,16),
                      gender = factor("Female", levels = c("Male", "Female")))
y_rep <- posterior_predict(womensrole_bglm_2, newdata)
summary(apply(y_rep, 1, diff))

## ---- rstanarm-rhat-fit, results='hide', warning=TRUE--------------------
bad_rhat <- stan_glm(mpg ~ ., data = mtcars, iter = 20, 
                     chains = CHAINS, cores = CORES, seed = SEED)
good_rhat <- update(bad_rhat, iter = 500, 
                    chains = CHAINS, cores = CORES, seed = SEED)

## ---- rstasnarm-rhat-bad-------------------------------------------------
rhat <- summary(bad_rhat)[, "Rhat"]
rhat[rhat > 1.1]

## ---- rstasnarm-rhat-good------------------------------------------------
any(summary(good_rhat)[, "Rhat"] > 1.1)

