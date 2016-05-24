## ----setup, echo = FALSE, results = "hide", message = FALSE--------------
set.seed(290875)

sapply(c("mlt", "survival", "eha", "prodlim", "truncreg", "lattice", "gridExtra",
         "MASS", "nnet", "HSAUR3", "sandwich", "flexsurv", "grid", "latticeExtra", 
         "colorspace", "multcomp"), library, char = TRUE)

if (!file.exists("analysis/DVC.rda")) {
    download.file("https://zenodo.org/record/17179/files/DVC.tgz", "DVC.tgz")
    untar("DVC.tgz", file = "analysis/DVC.rda", compressed = "gzip")
}
load("analysis/DVC.rda")
dvc <- c(margin.table(obs, 2))

logLik.phreg <- function(object) {
    ret <- object$loglik[2]
    attr(ret, "df") <- length(coef(object))
    class(ret) <- "logLik"
    ret
}
vcov.phreg <- function(object) object$var

trellis.par.set(list(plot.symbol = list(col=1,pch=20, cex=0.7),
                     box.rectangle = list(col=1),
                     box.umbrella = list(lty=1, col=1),
                     strip.background = list(col = "white")))
ltheme <- canonical.theme(color = FALSE)     ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme)

knitr::opts_chunk$set(echo = TRUE, results = 'markup', error = FALSE,
                      warning = FALSE, message = FALSE,
                      tidy = FALSE, cache = FALSE, size = "small",
                      fig.width = 6, fig.height = 4, fig.align = "center",
                      out.width = NULL, ###'.6\\linewidth', 
                      out.height = NULL,
                      fig.scap = NA)
knitr::render_sweave()  # use Sweave environments
knitr::set_header(highlight = '')  # do not \usepackage{Sweave}
## R settings
options(prompt = "R> ", continue = "+  ", useFancyQuotes = FALSE)  # JSS style
options(width = 75)

## ----geyser-var, echo = TRUE---------------------------------------------
library("mlt")
var_d <- numeric_var("duration", support = c(1.0, 5.0), 
                     add = c(-1, 1), bounds = c(0, Inf))

## ----geyser-basis, echo = TRUE-------------------------------------------
B_d <- Bernstein_basis(var = var_d, order = 8, ui = "increasing")

## ----geyser-ctm, echo = TRUE---------------------------------------------
ctm_d <- ctm(response = B_d, todistr = "Normal")

## ----geyser-grid, echo = TRUE--------------------------------------------
str(nd_d <- mkgrid(ctm_d, 200))

## ----geyser-data, echo = TRUE--------------------------------------------
data("geyser", package = "TH.data")
head(geyser)

## ----geyser-fit, echo = TRUE---------------------------------------------
mlt_d <- mlt(ctm_d, data = geyser)
logLik(mlt_d)

## ----geyser-density, echo = TRUE-----------------------------------------
nd_d$d <- predict(mlt_d, newdata = nd_d, type = "density")

## ----geyser-plot, echo = FALSE-------------------------------------------
plot(d ~ duration, data = nd_d, type = "l", ylab = "Density", xlab = "Duration time")

## ----CHFLS-1-------------------------------------------------------------
data("CHFLS", package = "HSAUR3")
polr_CHFLS_1 <- polr(R_happy ~ 1, data = CHFLS)

## ----CHFLS-1-basefun-----------------------------------------------------
nl <- nlevels(CHFLS$R_happy)
b_happy <- as.basis(~ R_happy, data = CHFLS, remove_intercept = TRUE,
                    contrasts.arg = list(R_happy = function(n) 
                        contr.treatment(n, base = nl)),
                    ui = diff(diag(nl - 1)), ci = rep(0, nl - 2))

## ----CHFLS-1-basefun-basis-----------------------------------------------
b_happy <- as.basis(CHFLS$R_happy)

## ----CHFLS-1-ctm---------------------------------------------------------
ctm_CHFLS_1 <- ctm(response = b_happy, todist = "Logistic")

## ----CHFLS-1-mlt---------------------------------------------------------
mlt_CHFLS_1 <- mlt(model = ctm_CHFLS_1, data = CHFLS)

## ----CHFLS-1-cmpr--------------------------------------------------------
logLik(polr_CHFLS_1)
logLik(mlt_CHFLS_1)
cbind(polr = polr_CHFLS_1$zeta, mlt = coef(mlt_CHFLS_1))

## ----CHFLS-1-pred--------------------------------------------------------
predict(polr_CHFLS_1, newdata = data.frame(1), type = "prob")
c(predict(mlt_CHFLS_1, newdata = data.frame(1), type = "density", 
          q = mkgrid(b_happy)[[1]]))
xtabs(~ R_happy, data = CHFLS) / nrow(CHFLS)

## ----geyser-w------------------------------------------------------------
var_w <- numeric_var("waiting", support = c(40.0, 100), add = c(-5, 15), 
                     bounds = c(0, Inf))
c(sapply(nd_w <- mkgrid(var_w, 100), range))

## ----geyser-bernstein----------------------------------------------------
B_w <- Bernstein_basis(var_w, order = 8, ui = "increasing")

## ----geyser-w-ctm, echo = TRUE-------------------------------------------
ctm_w <- ctm(response = B_w, todistr = "Normal")

## ----geyser-w-fit, echo = TRUE-------------------------------------------
mlt_w <- mlt(ctm_w, data = geyser)

## ----geyser-w-distribution, echo = TRUE----------------------------------
nd_w$d <- predict(mlt_w, newdata = nd_w, type = "distribution")

## ----geyser-w-plot, echo = FALSE-----------------------------------------
layout(matrix(1:2, ncol = 2))
plot(ecdf(geyser$waiting), col = "grey", xlab = "Waiting times", ylab = "Distribution", 
     main = "", cex = .75)
lines(nd_w$waiting, nd_w$d)
B_w_40 <- Bernstein_basis(order = 40, var = var_w, ui = "incre")
ctm_w_40 <- ctm(B_w_40, todistr = "Normal")
mlt_w_40 <- mlt(ctm_w_40, data = geyser, maxit = 5000)
nd_w$d2 <- predict(mlt_w_40, q = nd_w$waiting, type = "distribution")
lines(nd_w$waiting, nd_w$d2, lty = 2)
legend("bottomright", lty = 1:2, legend = c("M = 8", "M = 40"), bty = "n")
plot(nd_w$waiting, predict(mlt_w, q = nd_w$waiting, type = "density"), type = "l",
     ylim = c(0, .04), xlab = "Waiting times", ylab = "Density")
lines(nd_w$waiting, predict(mlt_w_40, q = nd_w$waiting, type = "density"), lty = 2)
rug(geyser$waiting, col = rgb(.1, .1, .1, .1))

## ----dvc-----------------------------------------------------------------
var_dvc <- numeric_var("dvc", support = min(dvc):max(dvc))
B_dvc <- Bernstein_basis(var_dvc, order = 6, ui = "increasing")
dvc_mlt <- mlt(ctm(B_dvc), data = data.frame(dvc = dvc))

## ----dvc-plot, echo = FALSE----------------------------------------------
q <- support(var_dvc)[[1]]
p <- predict(dvc_mlt, newdata = data.frame(1), q = q,
             type = "distribution")
plot(ecdf(dvc), col = "grey", xlab = "Number of Roe Deer-Vehicle Collisions",
     ylab = "Distribution", main = "", cex = .75)
lines(q, p, col = "blue")
lines(q, ppois(q, lambda = mean(dvc)), col = "darkred")
legend(400, .3, pch = c(20, NA, NA), lty = c(NA, 1, 1), 
       legend = c("ECDF", "Transformation Model", "Poisson"), bty = "n", cex = .8,
       col = c("grey", "blue", "darkred"))

## ----CHFLS-2, cache = TRUE-----------------------------------------------
polr_CHFLS_2 <- polr(R_happy ~ R_age + R_income, data = CHFLS)

## ----CHFLS-2-base--------------------------------------------------------
b_R <- as.basis(~ R_age + R_income, data = CHFLS, remove_intercept = TRUE, 
                negative = TRUE)

## ----CHFLS-2-ctm---------------------------------------------------------
ctm_CHFLS_2 <- ctm(response = b_happy, shifting = b_R, 
                   todistr = "Logistic")
mlt_CHFLS_2 <- mlt(ctm_CHFLS_2, data = CHFLS, scale = TRUE)

## ----CHFLS-2-cmpr--------------------------------------------------------
logLik(polr_CHFLS_2)
logLik(mlt_CHFLS_2)
cbind(polr = c(polr_CHFLS_2$zeta, coef(polr_CHFLS_2)), 
      mlt = coef(mlt_CHFLS_2))

## ----GBSG2-1, echo = TRUE------------------------------------------------
data("GBSG2", package = "TH.data")
GBSG2y <- numeric_var("y", support = c(100.0, max(GBSG2$time)), 
                      bounds = c(0, Inf))
GBSG2$y <- with(GBSG2, Surv(time, cens))

## ----GBSG2-1-Cox---------------------------------------------------------
B_GBSG2y <- Bernstein_basis(var = GBSG2y, order = 10, ui = "increasing")
fm_GBSG2 <- Surv(time, cens) ~ horTh + age + menostat + tsize + tgrade +
                               pnodes + progrec + estrec
ctm_GBSG2_1 <- ctm(B_GBSG2y, shifting = fm_GBSG2[-2L], data = GBSG2,
                   todistr = "MinExtrVal")

## ----GBSG2-1-xbasis, eval = FALSE----------------------------------------
#  as.basis(fm_GBSG2[-2L], data = GBSG2, remove_intercept)

## ----GBSG2-1-mlt---------------------------------------------------------
mlt_GBSG2_1 <- mlt(ctm_GBSG2_1, data = GBSG2, maxit = 3000, scale = TRUE)

## ----GBSG2-1-coxph-------------------------------------------------------
coxph_GBSG2_1 <- coxph(fm_GBSG2, data = GBSG2)
cf <- coef(coxph_GBSG2_1)
cbind(coxph = cf, mlt = coef(mlt_GBSG2_1)[names(cf)])

## ----GBSG2-1-fss, cache = TRUE-------------------------------------------
fss_GBSG2_1 <- flexsurvspline(fm_GBSG2, data = GBSG2, scale = "hazard", 
                              k = 9, bknots = log(support(GBSG2y)$y))
logLik(fss_GBSG2_1)
logLik(mlt_GBSG2_1)
cbind(coxph = cf, mlt = coef(mlt_GBSG2_1)[names(cf)],
      fss = coef(fss_GBSG2_1)[names(cf)])

## ----GBSG2-1-fss-plot, echo = FALSE--------------------------------------
p1 <- summary(fss_GBSG2_1, newdata = GBSG2[1,], ci = FALSE)
p2 <- predict(mlt_GBSG2_1, newdata = GBSG2[1, all.vars(fm_GBSG2[-2L])], 
              q = p1[[1]]$time, type = "survivor")
plot(p1[[1]]$time, p1[[1]]$est, type = "l", lty = 1, xlab = "Survival Time (days)", 
     ylab = "Probability", ylim = c(0, 1))
lines(p1[[1]]$time, p2[,1], lty = 2)
legend("topright", lty = 1:2, legend = c("flexsurvspline", "mlt"), bty = "n")

## ----GBSG2-2-------------------------------------------------------------
ly <- log_basis(GBSG2y, ui = "increasing")
ctm_GBSG2_2 <- ctm(ly, shifting = fm_GBSG2[-2L], data = GBSG2, 
                   negative = TRUE, todistr = "MinExtrVal")
mlt_GBSG2_2 <- mlt(ctm_GBSG2_2, data = GBSG2, fixed = c("log(y)" = 1), 
                   scale = TRUE)

## ----GBSG2-2-exp---------------------------------------------------------
survreg_GBSG2_2 <- survreg(fm_GBSG2, data = GBSG2, dist = "exponential")
phreg_GBSG2_2 <- phreg(fm_GBSG2, data = GBSG2, dist = "weibull", shape = 1)
logLik(survreg_GBSG2_2)
logLik(phreg_GBSG2_2)
logLik(mlt_GBSG2_2)
cbind(survreg = coef(survreg_GBSG2_2)[names(cf)], 
      phreg = -coef(phreg_GBSG2_2)[names(cf)], 
      mlt = coef(mlt_GBSG2_2)[names(cf)])

## ----GBSG2-3-------------------------------------------------------------
ctm_GBSG2_3 <- ctm(log_basis(GBSG2y, ui = "increasing"), 
                   shifting = fm_GBSG2[-2L], data = GBSG2, 
                   negative = TRUE, todistr = "MinExtrVal")
mlt_GBSG2_3 <- mlt(ctm_GBSG2_3, data = GBSG2, scale = TRUE)
survreg_GBSG2_3 <- survreg(fm_GBSG2, data = GBSG2, dist = "weibull")
phreg_GBSG2_3 <- phreg(fm_GBSG2, data = GBSG2, dist = "weibull")
logLik(survreg_GBSG2_3)
logLik(phreg_GBSG2_3)
logLik(mlt_GBSG2_3)
cbind(survreg = coef(survreg_GBSG2_3)[names(cf)] / survreg_GBSG2_3$scale, 
      phreg = - coef(phreg_GBSG2_3)[names(cf)], 
      mlt = coef(mlt_GBSG2_3)[names(cf)])

## ----BostonHousing-lm----------------------------------------------------
data("BostonHousing2", package = "mlbench")
lm_BH <- lm(cmedv ~ crim + zn + indus + chas + nox + rm + age + dis + 
            rad + tax + ptratio + b + lstat, data = BostonHousing2)

## ----BostonHousing-mlt---------------------------------------------------
BostonHousing2$medvc <- with(BostonHousing2, Surv(cmedv, cmedv < 50))
var_m <- numeric_var("medvc", support = c(10.0, 40.0), bounds = c(0, Inf))
fm_BH <- medvc ~ crim + zn + indus + chas + nox + rm + age + 
                 dis + rad + tax + ptratio + b + lstat

## ----BostonHousing-ctm---------------------------------------------------
B_m <- Bernstein_basis(var_m, order = 6, ui = "increasing")
ctm_BH <- ctm(B_m, shift = fm_BH[-2L], data = BostonHousing2, 
              todistr = "Normal")
mlt_BH <- mlt(ctm_BH, data = BostonHousing2, scale = TRUE)

## ----BostonHousing-plot, echo = FALSE, results = "hide"------------------
q <- 3:52
m <- predict(lm_BH, data = BostonHousing2)
s <- summary(lm_BH)$sigma
d <- sapply(m, function(m) pnorm(q, mean = m, sd = s))
nd <- expand.grid(q = q, lp = m)
nd$d <- c(d)

pfun <- function(x, y, z, subscripts, at, ...) {
    panel.contourplot(x, y, z, subscripts, at = 1:9/10, ...)
    panel.xyplot(x = m, y = BostonHousing2$cmedv, pch = 20,   
                 col = rgb(.1, .1, .1, .1),  ...)   
}
p1 <- contourplot(d ~ lp + q, data = nd, panel = pfun, 
                  xlab = "Linear predictor", ylab = "Observed", 
                  main = "Normal Linear Model")

d <- predict(mlt_BH, newdata = BostonHousing2, q = q, type="distribution")
lp <- c(predict(mlt_BH, newdata = BostonHousing2, q = 0, terms = "bshift"))
nd <- expand.grid(q = q, lp = -lp)
nd$d <- c(d)
pfun <- function(x, y, z, subscripts, at, ...) {
    panel.contourplot(x, y, z, subscripts, at = 1:9/10, ...)
    panel.xyplot(x = -lp, y = BostonHousing2$cmedv, pch = 20, 
                 col = rgb(.1, .1, .1, .1), ...)   
}
p2 <- contourplot(d ~ lp + q, data = nd, panel = pfun, xlab = "Linear predictor", ylab = "Observed", main = "Linear Transformation Model")
grid.arrange(p1, p2, nrow = 1)

## ----PSID1976------------------------------------------------------------
data("PSID1976", package = "AER")
PSID1976$nwincome <- with(PSID1976, (fincome - hours * wage)/1000)
PSID1976$hours <- as.double(PSID1976$hours)
PSID1976_0 <- subset(PSID1976, participation == "yes")
fm_PSID1976 <- hours ~ nwincome + education + experience + 
                       I(experience^2) + age + youngkids + oldkids

## ----PSID1976-truncreg---------------------------------------------------
tr_PSID1976 <- truncreg(fm_PSID1976, data = PSID1976_0)

## ----PSID1976-mlt--------------------------------------------------------
PSID1976_0$hours <- R(PSID1976_0$hours, tleft = 0)
b_hours <- as.basis(~ hours, data = PSID1976, 
                    ui = matrix(c(0, 1), nr  = 1), ci = 0)
ctm_PSID1976_1 <- ctm(b_hours, shift = fm_PSID1976[-2L], 
                      data = PSID1976_0, todistr = "Normal") 
mlt_PSID1976_1 <- mlt(ctm_PSID1976_1, data = PSID1976_0, scale = TRUE)

## ----PSID1976-cmpr-------------------------------------------------------
logLik(tr_PSID1976)
logLik(mlt_PSID1976_1)
cf <- coef(mlt_PSID1976_1)
cbind(truncreg = coef(tr_PSID1976),
      mlt = c(-cf[-grep("hours", names(cf))], 1) / cf["hours"])

## ----PSID1976-mlt-ctm----------------------------------------------------
var_h <- numeric_var("hours", support = range(PSID1976_0$hours$exact),
                     bounds = c(0, Inf))
B_hours <- Bernstein_basis(var_h, order = 6, ui = "increasing")
ctm_PSID1976_2 <- ctm(B_hours, shift = fm_PSID1976[-2L], 
                      data = PSID1976_0, todistr = "Normal")
mlt_PSID1976_2 <- mlt(ctm_PSID1976_2, data = PSID1976_0, 
                      scale = TRUE)
logLik(mlt_PSID1976_2)
AIC(mlt_PSID1976_1)
AIC(mlt_PSID1976_2)

## ----CHFLS-3, cache = TRUE-----------------------------------------------
b_health <- as.basis(~ R_health - 1, data = CHFLS)
ctm_CHFLS_3 <- ctm(b_happy, interacting = b_health, todist = "Logistic")
mlt_CHFLS_3 <- mlt(ctm_CHFLS_3, data = CHFLS, scale = TRUE,
                   maxit = 5000, gtol = 1e-3)
logLik(mlt_CHFLS_3)
predict(mlt_CHFLS_3, newdata = mkgrid(mlt_CHFLS_3), type = "distribution")

## ----CHFLS-4, cache = TRUE-----------------------------------------------
ctm_CHFLS_4 <- ctm(b_happy, interacting = b_health, shifting = b_R, 
                   todist = "Logistic")
mlt_CHFLS_4 <- mlt(ctm_CHFLS_4, data = CHFLS, scale = TRUE, 
                   maxit = 5000)
coef(mlt_CHFLS_4)[c("R_age", "R_income")]

## ----GBSG2-4-------------------------------------------------------------
b_horTh <- as.basis(GBSG2$horTh)
ctm_GBSG2_4 <- ctm(B_GBSG2y, interacting = b_horTh, 
                   todistr = "MinExtrVal")
mlt_GBSG2_4 <- mlt(ctm_GBSG2_4, data = GBSG2)

## ----GBSG2-strata-plot, echo = FALSE, results = "hide"-------------------
nd <- expand.grid(s <- mkgrid(mlt_GBSG2_4, 100))
nd$mlt_S <- c(predict(mlt_GBSG2_4, newdata = s, type = "survivor"))
nd$KM_S <- unlist(predict(prodlim(Surv(time, cens) ~ horTh, data = GBSG2), 
     	             newdata = data.frame(horTh = s$horTh), times = s$y))
plot(nd$y, nd$mlt_S, ylim = c(0, 1), xlab = "Survival time (days)",
     ylab = "Probability", type = "n", las = 1)
with(subset(nd, horTh == "no"), lines(y, mlt_S, col = "grey", lty = 2))
with(subset(nd, horTh == "yes"), lines(y, mlt_S, lty = 2))
with(subset(nd, horTh == "no"), lines(y, KM_S, type = "s", col = "grey"))
with(subset(nd, horTh == "yes"), lines(y, KM_S, type = "s"))
legend(250, 0.4, lty = c(1, 1, 2, 2), col = c("black", "grey", "black", "grey"),
       legend = c("hormonal therapy, KM", "no hormonal therapy, KM", 
                  "hormonal therapy, MLT", "no hormonal therapy, MLT"), bty = "n", cex = .75)

## ----GBSG2-5-------------------------------------------------------------
ctm_GBSG2_5 <- ctm(B_GBSG2y, interacting = b_horTh, shifting = ~ age, 
                   data = GBSG2, todistr = "MinExtrVal")
mlt_GBSG2_5 <- mlt(ctm_GBSG2_5, data = GBSG2, scale = TRUE)

## ----GBSG2-5-coxph-------------------------------------------------------
coxph_GBSG2_5 <- coxph(Surv(time, cens) ~ age + strata(horTh), 
                       data = GBSG2)
cf <- coef(coxph_GBSG2_5)
cbind(coxph = cf, mlt = coef(mlt_GBSG2_5)[names(cf)])

## ----CHFLS-5, cache = TRUE-----------------------------------------------
contrasts(CHFLS$R_health) <- "contr.treatment"
b_health <- as.basis(~ R_health, data = CHFLS)
ctm_CHFLS_5 <- ctm(b_happy, interacting = b_health, todist = "Logistic")
mlt_CHFLS_5 <- mlt(ctm_CHFLS_5, data = CHFLS, scale = TRUE, 
                   maxit = 10000, gtol = 1e-3)
predict(mlt_CHFLS_5, newdata = mkgrid(mlt_CHFLS_5), type = "distribution")
logLik(mlt_CHFLS_5)

## ----CHFLS-6, cache = TRUE-----------------------------------------------
b_R <- as.basis(~ R_age + R_income, data = CHFLS, remove_intercept = TRUE)
ctm_CHFLS_6 <- ctm(b_happy, interacting = b_R, todist = "Logistic")  
mlt_CHFLS_6 <- mlt(ctm_CHFLS_6, data = CHFLS, scale = TRUE,
                   maxit = 5000, gtol = 1e-3)
logLik(mlt_CHFLS_6)

## ----CHFLS-7, cache = TRUE-----------------------------------------------
ctm_CHFLS_7 <- ctm(b_happy, interacting = c(h = b_health, R = b_R), 
    todist = "Logistic")  
mlt_CHFLS_7 <- mlt(ctm_CHFLS_7, data = CHFLS, scale = TRUE,
                   maxit = 10000, gtol = 1e-3)
logLik(mlt_CHFLS_7)

## ----CHFLS-AIC-----------------------------------------------------------
c("1" = AIC(mlt_CHFLS_1), "2" = AIC(mlt_CHFLS_2), "3" = AIC(mlt_CHFLS_3),
  "4" = AIC(mlt_CHFLS_4), "5" = AIC(mlt_CHFLS_5), "6" = AIC(mlt_CHFLS_6),
  "7" = AIC(mlt_CHFLS_7))

## ----iris-1--------------------------------------------------------------
fm_iris <- Species ~ Sepal.Length + Sepal.Width + 
                     Petal.Length + Petal.Width
multinom_iris <- multinom(fm_iris, data = iris, trace = FALSE)
logLik(multinom_iris)
iris$oSpecies <- ordered(iris$Species)
b_Species <- as.basis(iris$oSpecies)
ctm_iris <- ctm(b_Species, 
                interacting = as.basis(fm_iris[-2L], data = iris), 
                todistr = "Logistic")
mlt_iris <- mlt(ctm_iris, data = iris, scale = TRUE)
logLik(mlt_iris)
p1 <- predict(mlt_iris, newdata = iris, q = sort(unique(iris$oSpecies)), 
              type = "density")
p2 <- predict(multinom_iris, newdata = iris, type = "prob")
max(abs(t(p1) - p2))

## ----GBSG2-6-------------------------------------------------------------
ctm_GBSG2_6 <- ctm(B_GBSG2y, shifting = ~ horTh, data = GBSG2, 
                   todistr = "MinExtrVal")
mlt_GBSG2_6 <- mlt(ctm_GBSG2_6, data = GBSG2)
logLik(mlt_GBSG2_6)
AIC(mlt_GBSG2_6)

## ----GBSG2-7-------------------------------------------------------------
b_horTh <- as.basis(~ horTh, data = GBSG2)
ctm_GBSG2_7 <- ctm(B_GBSG2y, interacting = b_horTh, 
                   todistr = "MinExtrVal")
mlt_GBSG2_7 <- mlt(ctm_GBSG2_7, data = GBSG2)
logLik(mlt_GBSG2_7)
AIC(mlt_GBSG2_7)

## ----GBSG2-deviation-plot, echo = FALSE, results = "hide"----------------

s <- mkgrid(mlt_GBSG2_7, 15)
s$y <- s$y[s$y > 100 & s$y < 2400]
nd <- expand.grid(s)
K <- model.matrix(ctm_GBSG2_7, data = nd)
Kyes <- K[nd$horTh == "yes",]
Kyes[,grep("Intercept", colnames(K))] <- 0  
gh <- glht(parm(coef(mlt_GBSG2_7), vcov(mlt_GBSG2_7)), Kyes)
ci <- confint(gh)
coxy <- s$y

K <- matrix(0, nrow = 1, ncol = length(coef(mlt_GBSG2_6)))
K[,length(coef(mlt_GBSG2_6))] <- 1
ci2 <- confint(glht(mlt_GBSG2_6, K))

plot(coxy, ci$confint[, "Estimate"], ylim = range(ci$confint), type = "n",
     xlab = "Survival time (days)", ylab = "Transformation deviation", las = 1)
polygon(c(coxy, rev(coxy)), c(ci$confint[,"lwr"], rev(ci$confint[, "upr"])),
        border = NA, col = rgb(.1, .1, .1, .1))
lines(coxy, ci$confint[, "Estimate"], lty = 1, lwd = 1)
lines(coxy, rep(ci2$confint[,"Estimate"], length(coxy)), lty = 2, lwd = 1) 
lines(coxy, rep(0, length(coxy)), lty = 3)
polygon(c(coxy[c(1, length(coxy))], rev(coxy[c(1, length(coxy))])),
        rep(ci2$confint[,c("lwr", "upr")], c(2, 2)),
        border = NA, col = rgb(.1, .1, .1, .1))
legend("bottomright", lty = 1:2, lwd = 1, legend = c("time-varying treatment effect",
       "time-constant log-hazard ratio"), bty = "n", cex = .75)

## ----GBSG2-8, echo = TRUE, cache = TRUE----------------------------------
var_a <- numeric_var("age", support = range(GBSG2$age))
B_age <- Bernstein_basis(var_a, order = 3)
b_horTh <- as.basis(GBSG2$horTh)
ctm_GBSG2_8 <- ctm(B_GBSG2y, 
                   interacting = b(horTh = b_horTh, age = B_age), 
                   todistr = "MinExtrVal")
mlt_GBSG2_8  <- mlt(ctm_GBSG2_8, data = GBSG2)
logLik(mlt_GBSG2_8)
AIC(mlt_GBSG2_8)

## ----GBSG2-8-plot, echo = FALSE------------------------------------------
nlev <- c(no = "without hormonal therapy", yes = "with hormonal therapy")
levels(nd$horTh) <- nlev[match(levels(nd$horTh), names(nlev))]
s <- mkgrid(mlt_GBSG2_8, 100)
nd <- expand.grid(s)
nd$s <- c(predict(mlt_GBSG2_8, newdata = s, type = "survivor"))
contourplot(s ~ age + y | horTh, data = nd, at = 1:9 / 10,
            ylab = "Survival time (days)", xlab = "Age (years)",
            scales = list(x = list(alternating = c(1, 1))))

## ----head, echo = TRUE, cache = TRUE-------------------------------------
data("db", package = "gamlss.data")
db$lage <- with(db, age^(1/3))
var_head = numeric_var("head", support = quantile(db$head, c(.1, .9)),
                       bounds = range(db$head))
B_head <- Bernstein_basis(var_head, order = 3, ui = "increasing")
var_lage <- numeric_var("lage", support = quantile(db$lage, c(.1, .9)),
                        bounds = range(db$lage))
B_age <- Bernstein_basis(var_lage, order = 3, ui = "none")
ctm_head <- ctm(B_head, interacting = B_age)
mlt_head <- mlt(ctm_head, data = db, maxit = 5000, scale = TRUE)

## ----head-plot, echo = FALSE---------------------------------------------
pr <- expand.grid(s <- mkgrid(ctm_head, 100))
pr$p <- c(predict(mlt_head, newdata = s, type = "distribution"))
pr$lage <- pr$lage^3
pr$cut <- factor(pr$lage > 2.5)
levels(pr$cut) <- c("Age < 2.5 yrs", "Age > 2.5 yrs")
pfun <- function(x, y, z, subscripts, at, ...) {
    panel.contourplot(x, y, z, subscripts,
        at = c(0.4, 2, 10, 25, 50, 75, 90, 98, 99.6)/ 100, ...)
    panel.xyplot(x = db$age, y = db$head, pch = 20,
                 col = rgb(.1, .1, .1, .1), ...)
}
print(contourplot(p ~ lage + head | cut, data = pr, panel = pfun, region = FALSE,
            xlab = "Age (years)", ylab = "Head circumference (cm)",
            scales = list(x = list(relation = "free"))))

## ----BostonHousing-dr, cache = TRUE--------------------------------------
b_BH_s <- as.basis(fm_BH[-2L], data = BostonHousing2, scale = TRUE)
ctm_BHi <- ctm(B_m, interacting = b_BH_s, sumconstr = FALSE)
mlt_BHi <- mlt(ctm_BHi, data = BostonHousing2, scale = TRUE, 
               maxit = 5000, gtol = 1e-3)
logLik(mlt_BHi)
AIC(mlt_BHi)
AIC(mlt_BH)

## ----Boston-Housing-dr-plot, echo = FALSE--------------------------------
q <- mkgrid(var_m, 100)[[1]]
tr <- predict(mlt_BH, newdata = BostonHousing2[, all.vars(fm_BH[-2L])],
              q = q, type = "density")
tri <- predict(mlt_BHi, newdata = BostonHousing2[, all.vars(fm_BH[-2L])],
              q = q, type = "density")
layout(matrix(1:2, ncol = 2))
Q <- matrix(q, nrow = length(q), ncol = ncol(tr))
ylim <- range(c(tr, tri))
matplot(Q, tr, ylim = ylim, xlab = "Median Value", ylab = "Density",
        type = "l", col = rgb(.1, .1, .1, .1), lty = 1) 
matplot(Q, tri, ylim = ylim, xlab = "Median Value", ylab = "Density",
        type = "l", col = rgb(.1, .1, .1, .1), lty = 1)

## ----treepipit, echo = TRUE, cache = TRUE--------------------------------
data("treepipit", package = "coin")
treepipit$ocounts <- ordered(treepipit$counts)
B_cs <- Bernstein_basis(var = numeric_var("coverstorey", support = 1:110), 
                        order = 4)
B_c <- as.basis(treepipit$ocounts)
ctm_treepipit <- ctm(B_c, interacting = B_cs)
mlt_treepipit <- mlt(ctm_treepipit, data = treepipit, maxit = 10000, 
                     gtol = 1e-3)

## ----mgcv, echo = FALSE, results = "hide"--------------------------------
library("mgcv") ### masks nnet::multinom

## ----treepipit-gam-------------------------------------------------------
gam_treepipit <- gam(counts ~ s(coverstorey), data = treepipit, 
                     family = "poisson")

## ----treepipit-plot, echo = FALSE, fig.height = 3------------------------
s <- mkgrid(ctm_treepipit, 100)
s$ocounts <- s$ocounts[1:5]
nd <- expand.grid(s)
nd$p <- c(predict(mlt_treepipit, newdata = s, type = "distribution"))

### produce a table
tpt <- xtabs(~ counts + coverstorey, data = treepipit)

### construct a data frame with frequencies
treepipit2 <- sapply(as.data.frame(tpt, stringsAsFactors = FALSE),
                     as.integer)

s <- mkgrid(ctm_treepipit, 10)
s$ocounts <- s$ocounts[1]
K <- model.matrix(ctm_treepipit, data = expand.grid(s))
#g <- glht(parm(coef(mod), vcov(mod)), linfct = K)
#confint(g)
nd$lambda <- predict(gam_treepipit, newdata = nd, type = "response")

layout(matrix(1:3, nr = 1))
par("mai" = par("mai") * c(1, .95, 1, .85))
xlim <- range(treepipit[, "coverstorey"]) * c(0.98, 1.05)
xlab <- "Cover storey"
ylab <- "Number of tree pipits (TP)"
### scatterplot again; plots are proportional to number of plots
plot(counts ~ coverstorey, data = treepipit2, cex = sqrt(Freq),
     ylim = c(-.5, 5), xlab = xlab, ylab = ylab, col = "darkgrey", 
     xlim = xlim, las = 1, main = "Observations")

plot(c(0, 110), c(0, 1), type = "n", xlab = xlab, ylab = "Conditional probability",
     xlim = xlim, las = 1, main = "MLT")
with(subset(nd, ocounts == "0"), lines(coverstorey, p, lty = 1))
with(subset(nd, ocounts == "1"), lines(coverstorey, p, lty = 2))
with(subset(nd, ocounts == "2"), lines(coverstorey, p, lty = 3))
with(subset(nd, ocounts == "3"), lines(coverstorey, p, lty = 4))
with(subset(nd, ocounts == "4"), lines(coverstorey, p, lty = 5))
abline(h = 1, lty = 6)
legend("bottomright", lty = 1:6, legend = c(expression(TP == 0),
                                            expression(TP <= 1),
                                            expression(TP <= 2),
                                            expression(TP <= 3),
                                            expression(TP <= 4),
                                            expression(TP <= 5)), bty = "n")

plot(c(0, 110), c(0, 1), type = "n", xlab = xlab, ylab = "Conditional probability",
     xlim = xlim, las = 1, main = "GAM")
with(subset(nd, ocounts == "0"), lines(coverstorey, ppois(0, lambda), lty = 1))
with(subset(nd, ocounts == "1"), lines(coverstorey, ppois(1, lambda), lty = 2))
with(subset(nd, ocounts == "2"), lines(coverstorey, ppois(2, lambda), lty = 3))
with(subset(nd, ocounts == "3"), lines(coverstorey, ppois(3, lambda), lty = 4))
with(subset(nd, ocounts == "4"), lines(coverstorey, ppois(4, lambda), lty = 5))
abline(h = 1, lty = 6)

## ----CHFLS-2-cmpr-2------------------------------------------------------
max(abs(estfun(polr_CHFLS_2) - (-estfun(mlt_CHFLS_2)[,c(4, 5, 1:3)])))
cbind(polr = sqrt(diag(vcov(polr_CHFLS_2))),
      mlt = sqrt(diag(vcov(mlt_CHFLS_2)))[c(4, 5, 1:3)])
cftest(polr_CHFLS_2)
cftest(mlt_CHFLS_2, parm = names(coef(polr_CHFLS_2)))

## ----GBSG2-1-coxph-cmpr--------------------------------------------------
cf <- coef(coxph_GBSG2_1)
cbind(coxph = sqrt(diag(vcov(coxph_GBSG2_1))),
      mlt = sqrt(diag(vcov(mlt_GBSG2_1)))[names(cf)],
      fss = sqrt(diag(vcov(fss_GBSG2_1)))[names(cf)])

## ----GBSG2-1-coxph-cmpr-cftest-------------------------------------------
cftest(coxph_GBSG2_1)
cftest(mlt_GBSG2_1, parm = names(cf))
cftest(fss_GBSG2_1, parm = names(cf))

## ----geyser-w-band-------------------------------------------------------
cb_w <- confband(mlt_w, newdata = data.frame(1), K = 20, cheat = 100)

## ----geyser-w-cbplot, echo = FALSE---------------------------------------
layout(matrix(1:2, ncol = 2))
#i <- (cb_w[, "q"] > 45 & cb_w[, "q"] < 110)
#cb_w[-i, "lwr"] <- NA
#cb_w[-i, "upr"] <- NA
plot(cb_w[, "q"], cb_w[, "Estimate"], xlab = "Waiting times", ylab = "Transformation", 
     main = "", type = "l")

q <- cb_w[, "q"]
lwr <- cb_w[, "lwr"]
upr <-  cb_w[, "upr"]
polygon(c(q, rev(q)), c(lwr, rev(upr)),
        border = NA, col = rgb(.1, .1, .1, .1))

#lines(cb_w[, "q"], cb_w[, "lwr"])
#lines(cb_w[, "q"], cb_w[, "upr"])
rug(geyser$waiting, col = rgb(.1, .1, .1, .1))
plot(ecdf(geyser$waiting), col = "grey", xlab = "Waiting times", ylab = "Distribution", 
     main = "", cex = .5)
lines(cb_w[, "q"], pnorm(cb_w[, "Estimate"]))

polygon(c(q, rev(q)), c(pnorm(lwr), rev(pnorm(upr))),                                     
        border = NA, col = rgb(.1, .1, .1, .1))

# lines(cb_w[, "q"], pnorm(cb_w[, "lwr"]))
# lines(cb_w[, "q"], pnorm(cb_w[, "upr"]))
rug(geyser$waiting, col = rgb(.1, .1, .1, .1))

## ----geyser-w-simulate, results = "hide", cache = TRUE-------------------
new_w <- simulate(mlt_w, nsim = 100)
llr <- numeric(length(new_w))
pdist <- vector(mode = "list", length = length(new_w))
pdens <- vector(mode = "list", length = length(new_w))
ngeyser <- geyser
q <- mkgrid(var_w, 100)[[1]]
for (i in 1:length(new_w)) {
    ngeyser$waiting <- new_w[[i]]
    mlt_i <- mlt(ctm_w, data = ngeyser, scale = TRUE, 
                 theta = coef(mlt_w))
    llr[[i]] <- logLik(mlt_i) - logLik(mlt_i, parm = coef(mlt_w))
    pdist[[i]] <- predict(mlt_i, newdata = data.frame(1), 
                          type = "distribution", q = q)
    pdens[[i]] <- predict(mlt_i, newdata = data.frame(1), 
                          type = "density", q = q)
}

## ----geyser-w-simulate-plot, echo = FALSE--------------------------------
i <- which(llr < quantile(llr, prob = .95))
tpdist <- pdist[i]
tpdens <- pdens[i]
layout(matrix(1:2, ncol = 2))
plot(q, tpdist[[1]], type = "n", xlab = "Waiting times", ylab = "Distribution")
polygon(c(q, rev(q)), c(pnorm(lwr), rev(pnorm(upr))),
        border = NA, col = "greenyellow")
tmp <- sapply(tpdist, function(x) lines(q, x, col = rgb(.1, .1, .1, .1)))
plot(q, tpdens[[1]], type = "n", ylim = range(unlist(pdens)),
     xlab = "Waiting times", ylab = "Density")
tmp <- sapply(tpdens, function(x) lines(q, x, col = rgb(.1, .1, .1, .1)))

## ----variables-factor----------------------------------------------------
f_eye <- factor_var("eye", desc = "eye color", 
                    levels = c("blue", "brown", "green", "grey", "mixed"))

## ----variables-factor-methods--------------------------------------------
variable.names(f_eye)
desc(f_eye)
variables::unit(f_eye)
support(f_eye)
bounds(f_eye)
is.bounded(f_eye)

## ----variables-factor-mkgrid---------------------------------------------
mkgrid(f_eye)

## ----variables-ordered---------------------------------------------------
o_temp <- ordered_var("temp", desc = "temperature", 
                      levels = c("cold", "lukewarm", "warm", "hot"))

## ----variables-ordered-methods-------------------------------------------
variable.names(o_temp)
desc(o_temp)
variables::unit(o_temp)
support(o_temp)
bounds(o_temp) 
is.bounded(o_temp)
mkgrid(o_temp)

## ----variables-fd--------------------------------------------------------
v_age <- numeric_var("age", desc = "age of patient", 
                     unit = "years", support = 25:75)

## ----variables-fd-methods------------------------------------------------
variable.names(v_age)
desc(v_age)
variables::unit(v_age)
support(v_age)
bounds(v_age) 
is.bounded(v_age)

## ----variables-fd-mkgrid-------------------------------------------------
mkgrid(v_age)

## ----variables-c---------------------------------------------------------
v_temp <- numeric_var("ztemp", desc = "Zurich daytime temperature", 
                      unit = "Celsius", support = c(-10.0, 35.0), 
                      add = c(-5, 5), bounds = c(-273.15, Inf))

## ----variables-c-methods-------------------------------------------------
variable.names(v_temp)
desc(v_temp)
variables::unit(v_temp)
support(v_temp)
bounds(v_temp) 
is.bounded(v_temp)

## ----variables-c-mkgrid--------------------------------------------------
mkgrid(v_temp, n = 20)

## ----variables-vars------------------------------------------------------
vars <- c(f_eye, o_temp, v_age, v_temp)

## ----variables-vars-methods----------------------------------------------
variable.names(vars)
desc(vars) 
variables::unit(vars)
support(vars)
bounds(vars)
is.bounded(vars)
mkgrid(vars, n = 20)

## ----variables-vars-expand-----------------------------------------------
nd <- expand.grid(mkgrid(vars))

## ----variables-check-----------------------------------------------------
check(vars, data = nd)

## ----basefun-polynom-----------------------------------------------------
xvar <- numeric_var("x", support = c(0.1, pi), bounds= c(0, Inf))
x <- as.data.frame(mkgrid(xvar, n = 20))
### set-up basis of order 3 ommiting the quadratic term
class(pb <- polynomial_basis(xvar, coef = c(TRUE, TRUE, FALSE, TRUE)))

## ----basefun-polynom-fun-------------------------------------------------
head(pb(x))

## ----basefun-polynom-mm--------------------------------------------------
head(model.matrix(pb, data = x))

## ----basefun-polynom-pred------------------------------------------------
### evaluate polynomial defined by basis and coefficients
predict(pb, newdata = x, coef = c(1, 2, 0, 1.75))
### evaluate 1st derivative
predict(pb, newdata = x, coef = c(1, 2, 0, 1.75), deriv = c(x = 1L))

## ----basefun-log---------------------------------------------------------
### set-up log-basis with intercept for positive variable
class(lb <- log_basis(xvar, ui = "increasing"))
head(X <- model.matrix(lb, data = x))

## ----basefun-log-constr--------------------------------------------------
attr(X, "constraint")

## ----basefun-log-pred----------------------------------------------------
predict(lb, newdata = x, coef = c(1, 2))
predict(lb, newdata = x, coef = c(1, 2), deriv = c(x = 1L))

## ----basefun-Bernstein---------------------------------------------------
class(bb <- Bernstein_basis(xvar, order = 3, ui = "increasing"))
head(X <- model.matrix(bb, data = x))

## ----basefun-Bernstein-constr--------------------------------------------
cf <- c(1, 2, 2.5, 2.6)
(cnstr <- attr(X, "constraint"))
all(cnstr$ui %*% cf > cnstr$ci)

## ----basefun-Bernstein-predict-------------------------------------------
predict(bb, newdata = x, coef = cf)
predict(bb, newdata = x, coef = cf, deriv = c(x = 1))

## ----basefun-as.basis----------------------------------------------------
iv <- as.vars(iris)
fb <- as.basis(~ Species + Sepal.Length + Sepal.Width,  data = iv,
               remove_intercept = TRUE, negative = TRUE, 
               contrasts.args =  list(Species = "contr.sum"))
class(fb)
head(model.matrix(fb, data = iris))

## ----basefun-c-----------------------------------------------------------
class(blb <- c(bern = bb, 
               log = log_basis(xvar, ui = "increasing", 
                               remove_intercept = TRUE)))
head(X <- model.matrix(blb, data = x))
attr(X, "constraint")

## ----basefun-b-----------------------------------------------------------
fb <- as.basis(~ g, data = factor_var("g", levels = LETTERS[1:2]))
class(bfb <- b(bern = bb, f = fb))
nd <- expand.grid(mkgrid(bfb, n = 10))
head(X <- model.matrix(bfb, data = nd))
attr(X, "constraint")

## ----basefun-b-sumconstr-------------------------------------------------
bfb <- b(bern = bb, f = fb, sumconstr = TRUE)
head(X <- model.matrix(bfb, data = nd))
attr(X, "constraint")

## ----R-ordered-----------------------------------------------------------
head(R(sort(unique(CHFLS$R_happy))))

## ----R-integer-----------------------------------------------------------
R(1:5, bounds = c(1, 5))

## ----R-numeric-----------------------------------------------------------
x <- rnorm(10)
xt <- round(x[x > -1 & x <= 2], 3)
xl <- xt - sample(c(Inf, (0:(length(xt) - 2)) / length(xt)), 
                  replace = FALSE)
xr <- xt + sample(c(Inf, (0:(length(xt) - 2)) / length(xt)), 
                  replace = FALSE)
R(c(1.2, rep(NA, length(xt))), cleft = c(NA, xl), cright = c(NA, xr), 
  tleft = -1, tright = 2)

## ----R-Surv--------------------------------------------------------------
head(geyser$duration)
head(R(geyser$duration))

## ----mlt-chi-p-----------------------------------------------------------
pY <- function(x) pchisq(x, df = 20)
dY <- function(x) dchisq(x, df = 20)
qY <- function(p) qchisq(p, df = 20)

## ----mlt-chi-B-----------------------------------------------------------
yvar <- numeric_var("y", support = qY(c(.001, 1 - .001)), 
                    bounds = c(0, Inf))
By <- Bernstein_basis(yvar, order = ord <- 15, ui = "increasing")

## ----mlt-chi-mlt---------------------------------------------------------
m <- ctm(By)
d <- as.data.frame(mkgrid(yvar, n = 500))
mod <- mlt(m, data = d, dofit = FALSE)

## ----mlt-chi-trafo-------------------------------------------------------
h <- function(x) qnorm(pY(x))
x <- seq(from = support(yvar)[["y"]][1], to = support(yvar)[["y"]][2], 
         length.out = ord + 1)
mlt::coef(mod) <- h(x)

## ----mlt-chi-sim---------------------------------------------------------
d$grid <- d$y
d$y <- simulate(mod)
fmod <- mlt(m, data = d, scale = TRUE)

## ----mlt-chi-model-------------------------------------------------------
coef(mod)
coef(fmod)
logLik(fmod)
logLik(fmod, parm = coef(mod))

## ----mlt-chi-plot--------------------------------------------------------
## compute true density
d$dtrue <- dY(d$grid)
d$dest <- predict(fmod, q = sort(d$grid), type = "density")
plot(mod, newdata = d, type = "density", col = "black", 
     xlab = "y", ylab = "Density", ylim = c(0, max(d$dest)))
lines(d$grid, d$dtrue, lty = 2)
lines(sort(d$grid), d$dest[order(d$grid)], lty = 3)
legend("topright", lty = 1:3, bty = "n", 
       legend = c("True", "Approximated", "Estimated"))

## ----mlt-coef, echo = FALSE, results = "hide"----------------------------
### print coefs for regression tests
objs <- ls()
mltobj <- objs[grep("^mlt_", objs)]
sapply(mltobj, function(m) eval(parse(text = paste("coef(", m, ")"))))
#library("mlt.docreg")
#sapply(mltobj, function(m) eval(parse(text = paste("checkGH(", m, ")"))))

## ----sessionInfo, echo = FALSE, results = "hide"-------------------------
sessionInfo()

