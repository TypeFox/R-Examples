## ----include=FALSE-------------------------------------------------------
library("knitr")
options(prompt = "R> ", continue = "+  ", width = 75, useFancyQuotes = FALSE)
opts_chunk$set(fig.align = 'center')
render_sweave()
library("KFAS")

## ----'gaussian_model'----------------------------------------------------
data("alcohol")
deaths <- window(alcohol[, 2], end = 2007)
population <- window(alcohol[, 6], end = 2007)
Zt <- matrix(c(1, 0), 1, 2)
Ht <- matrix(NA)
Tt <- matrix(c(1, 0, 1, 1), 2, 2)
Rt <- matrix(c(1, 0), 2, 1)
Qt <- matrix(NA)
a1 <- matrix(c(1, 0), 2, 1)
P1 <- matrix(0, 2, 2)
P1inf <- diag(2)

model_gaussian <- SSModel(deaths / population ~ -1 +
    SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1,
      P1 = P1, P1inf = P1inf), H = Ht)

## ----'gaussian_fit'------------------------------------------------------
fit_gaussian <- fitSSM(model_gaussian, inits = c(0, 0), method = "BFGS")
out_gaussian <- KFS(fit_gaussian$model)

## ----'gaussian_plot', fig.pos = '!ht', fig.cap = 'Alcohol-related deaths in Finland in the age group of 40--49 years (black line) with predicted (red) and smoothed (blue) estimates.',out.width='0.7\\linewidth', echo = FALSE----
ts.plot(cbind(deaths / population,
  out_gaussian$a[-(length(deaths) + 1), 1], out_gaussian$alphahat[, 1]),
  ylab = "Alcohol-related deaths in Finland per 100,000 persons",
  xlab = "Year", col = c(1, 2, 4))

## ----'nongaussian_model'-------------------------------------------------
model_poisson <- SSModel(deaths ~ -1 +
    SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, P1inf = P1inf),
  distribution = "poisson", u = population)

## ----'nongaussian_example_fit', echo = FALSE-----------------------------
fit_poisson <- fitSSM(model_poisson, inits = -5, method = "BFGS")
out_poisson <- KFS(fit_poisson$model, smoothing = "state")

## ----'nongaussian_example', fig.pos = '!ht', fig.cap = 'Alcohol-related deaths in Finland (black line) with smoothed estimates from Gaussian model (blue) and Poisson model (red).',out.width='0.7\\linewidth', echo = FALSE----
ts.plot(deaths / population, out_gaussian$alphahat[, 1], exp(out_poisson$alphahat[, 1]),
  ylab = "Alcohol-related deaths in Finland per 100,000 persons",
  xlab = "Year", col = c(1,4,2))

## ----'structural_time_series'--------------------------------------------
model_structural <- SSModel(deaths / population ~
    SSMtrend(degree = 2, Q = list(matrix(NA), matrix(0))), H = matrix(NA))
fit_structural <- fitSSM(model_structural,
  inits = c(0, 0), method = "BFGS")
fit_structural$model["Q"]

## ----'arima_time_series'-------------------------------------------------
drift <- 1:length(deaths)
model_arima <- SSModel(deaths / population ~ drift +
    SSMarima(ma = 0, d = 1, Q = 1))

update_model <- function(pars, model) {
  tmp <- SSMarima(ma = pars[1], d = 1, Q = pars[2])
  model["R", states = "arima"] <- tmp$R
  model["Q", states = "arima"] <- tmp$Q
  model["P1", states = "arima"] <- tmp$P1
  model
}

fit_arima <- fitSSM(model_arima, inits = c(0, 1),
  updatefn = update_model, method = "L-BFGS-B",
  lower = c(-1, 0), upper = c(1, 100))
fit_arima$optim.out$par

## ----'arima_structural_comparison'---------------------------------------
(out_arima <- KFS(fit_arima$model))
(out_structural <- KFS(fit_structural$model))
out_arima$logLik
out_structural$logLik

## ----'glmexample1'-------------------------------------------------------
counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
model_glm1 <- SSModel(counts ~ outcome + treatment,
  distribution = "poisson")

## ----'glmexample2'-------------------------------------------------------
model_glm2 <- SSModel(counts ~ SSMregression(~ outcome + treatment),
  distribution = "poisson")

## ----'lmmexample1'-------------------------------------------------------
library("lme4", quietly = TRUE)
y_split <- split(sleepstudy["Reaction"], sleepstudy["Subject"])
p <- length(y_split)
y <- matrix(unlist(y_split), ncol = p,
  dimnames = list(NULL, paste("Subject", names(y_split))))

## ----'lmmexample2'-------------------------------------------------------
dataf <- split(sleepstudy, sleepstudy["Subject"])

## ----'lmmexample4'-------------------------------------------------------
P1 <- as.matrix(.bdiag(replicate(p, matrix(NA, 2, 2), simplify = FALSE)))
model_lmm <- SSModel(y ~ - 1 +
    SSMregression(rep(list(~ Days), p), type = "common", data = dataf,
      remove.intercept = FALSE) +
    SSMregression(rep(list(~ Days), p), data = dataf,
      remove.intercept = FALSE, P1 = P1),
  H = diag(NA, p))

## ----'lmmexample5'-------------------------------------------------------
model_lmm2 <- SSModel(y ~ - 1 +
    SSMregression(~ Days, type = "common", remove.intercept = FALSE) +
    SSMregression(~ Days, remove.intercept = FALSE, P1 = P1),
  H = diag(NA, p), data = data.frame(Days = 0:9))

## ----'estimate_lmm'------------------------------------------------------
update_lmm <- function(pars, model) {
  P1 <- diag(exp(pars[1:2]))
  P1[1, 2] <- pars[3]
  P1 <- crossprod(P1)
  model["P1", states = 3:38] <-
    as.matrix(.bdiag(replicate(p, P1, simplify = FALSE)))
  model["H"] <- diag(exp(pars[4]), p)
  model
}

fit_lmm <- fitSSM(model_lmm, c(1, 1, 1, 5), update_lmm, method = "BFGS")

## ----'customexample_model'-----------------------------------------------
model_poisson <- SSModel(deaths ~ SSMtrend(2, Q = list(NA, 0)) +
    SSMcustom(Z = 1, T = 0, Q = NA, P1 = NA),
  distribution = "poisson",  u = population)

## ----'customexample_fit'-------------------------------------------------
update_poisson <- function(pars, model) {
  model["Q", etas = "level"] <- exp(pars[1])
  model["Q", etas = "custom"] <- exp(pars[2])
  model["P1", states = "custom"] <- exp(pars[2])
  model
}
fit_poisson <- fitSSM(model_poisson, c(-3, -3),
  update_poisson, method = "BFGS")
fit_poisson$model["Q", etas = "level"]
fit_poisson$model["Q", etas = "custom"]

## ----'customexample_output',fig.pos = '!ht', fig.cap = 'Alcohol-related deaths in Finland (black line) with smoothed estimates from Gaussian model (blue) and Poisson model with additional noise (red).',out.width='0.7\\linewidth', echo =FALSE----
out_poisson <- KFS(fit_poisson$model)
ts.plot(deaths / population, coef(out_structural, states = "level"),
  exp(coef(out_poisson, states = "level")),
  ylab = "Alcohol-related deaths in Finland per 100,000 persons",
  xlab = "Year", col = c(1, 4, 2))

## ----'alcoholPlot1', fig.pos = '!ht', fig.cap = 'Alcohol-related deaths per 100,000 persons in Finland in 1969--2007 for four age groups.',out.width='0.7\\linewidth'----
data("alcohol")
colnames(alcohol)
ts.plot(window(alcohol[, 1:4] / alcohol[, 5:8], end = 2007), col = 1:4,
  ylab = "Alcohol-related deaths in Finland per 100,000 persons",
  xlab = "Year")
legend("topleft",col = 1:4, lty = 1, legend = colnames(alcohol)[1:4])

## ----'alcoholfit1'-------------------------------------------------------
alcoholPred <- window(alcohol, start = 1969, end = 2007)
model <- SSModel(alcoholPred[, 1:4] ~
    SSMtrend(2, Q = list(matrix(NA, 4, 4), matrix(0, 4, 4))) +
    SSMcustom(Z = diag(1, 4), T = diag(0, 4), Q = matrix(NA, 4, 4),
      P1 = matrix(NA, 4, 4)), distribution = "poisson",
  u = alcoholPred[, 5:8])

## ----'alcoholfit2'-------------------------------------------------------
updatefn <- function(pars, model, ...){
  Q <- diag(exp(pars[1:4]))
  Q[upper.tri(Q)] <- pars[5:10]
  model["Q", etas = "level"] <- crossprod(Q)
  Q <- diag(exp(pars[11:14]))
  Q[upper.tri(Q)] <- pars[15:20]
  model["Q", etas = "custom"] <- model["P1", states = "custom"] <-
    crossprod(Q)
  model
}

## ----'alcoholfit3'-------------------------------------------------------
init <- chol(cov(log(alcoholPred[, 1:4] / alcoholPred[, 5:8])) / 10)
fitinit <- fitSSM(model, updatefn = updatefn,
  inits = rep(c(log(diag(init)), init[upper.tri(init)]), 2),
  method = "BFGS")
-fitinit$optim.out$val

fit <- fitSSM(model, updatefn = updatefn, inits = fitinit$optim.out$par,
  method = "BFGS", nsim = 250)
-fit$optim.out$val

## ----'alcoholfit4'-------------------------------------------------------
varcor <- fit$model["Q", etas = "level"]
varcor[upper.tri(varcor)] <- cov2cor(varcor)[upper.tri(varcor)]
print(varcor, digits = 2)
varcor <- fit$model["Q", etas = "custom"]
varcor[upper.tri(varcor)] <- cov2cor(varcor)[upper.tri(varcor)]
print(varcor, digits = 2)

## ----'KFS'---------------------------------------------------------------
out <- KFS(fit$model, nsim = 1000)
out

## ----'states',fig.pos = '!ht', fig.cap = 'Smoothed level and white noise components.',out.width='\\linewidth'----
plot(coef(out, states = c("level", "custom")), main = "Smoothed states",
  yax.flip = TRUE)

## ----'diagnostics1',fig.pos = '!ht', fig.cap = 'Autocorrelations and cross-correlations of recursive residuals.',out.width='\\linewidth'----
res <- rstandard(KFS(fit$model, filtering = "mean", smoothing = "none",
  nsim = 1000))
acf(res, na.action = na.pass)

## ----'prediction'--------------------------------------------------------
pred <- predict(fit$model,
  newdata = SSModel(ts(matrix(NA, 6, 4), start = 2008) ~ - 1 +
      SSMcustom(Z = fit$model$Z, T = fit$model$T, R = fit$model$R,
        Q = fit$model$Q), u = 1, distribution= "poisson"),
  interval = "confidence", nsim = 10000)

## ----'predictplot',fig.pos = '!ht', fig.cap = 'Observed number of alcohol related deaths per 100,000 persons in Finland (black), fitted values (red) and intensity predictions for years the 2008--2013 together with 95\\% prediction intervals (green).',out.width='\\linewidth'----
trend <- exp(signal(out, "trend")$signal)
par(mfrow = c(2, 2), mar = c(2, 2, 2, 2) + 0.1, oma = c(2, 2, 0, 0))
for(i in 1:4)
  ts.plot(alcohol[, i]/alcohol[, 4+i], trend[, i], pred[[i]],
    col = c(1, 2, rep(3, 3)), xlab = NULL, ylab = NULL,
    main = colnames(alcohol)[i])
mtext("Number of alcohol related deaths per 100,000 persons in Finland",
  side = 2, outer = TRUE)
mtext("Year", side = 1, outer = TRUE)

## ----'salmonellaINLA'----------------------------------------------------
#   library("INLA", quietly = TRUE)
# data("Salm")
# mod.salm <- inla(y ~ log(dose + 10) + dose +
#     f(rand, model = "iid", param = c(0.001, 0.001)),
#   family = "poisson", data = Salm)
# h.salm <- inla.hyperpar(mod.salm)

## ----'salmonellaKFAS'----------------------------------------------------
#   Salm$rand <- as.factor(Salm$rand)
# model <- SSModel(y ~ log(dose + 10) + dose +
#     SSMregression(~ - 1 + rand, P1 = diag(NA, 18),
#       remove.intercept = FALSE),
#   data = Salm, distribution = "poisson")
#
# updatefn<-function(pars,model,...){
#   diag(model["P1", states = 4:21]) <- exp(pars)
#   model
# }
#
# fit <- fitSSM(model, updatefn = updatefn, inits = -3, method = "BFGS",
#   nsim = 1000)

## ----'salmonellaResults'-------------------------------------------------
#   out <- KFS(fit$model, nsim = 10000)
# out
# h.salm$summary.fixed[, 1:2]
# h.salm$summary.random$rand[, 2:3]
# 1 / h.salm$summary.hyper[1]
# fit$model["P1", states = 4]

