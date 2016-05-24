## MuMIn tests: support for various model classes

require(MuMIn)
.checkPkg <- function(package) 
length(find.package(package, quiet = TRUE)) == length(package)


options(na.action = "na.fail")

# --------------------------------------------------------------------------------
# TEST gls
library(nlme) 

fm1Dial.gls <- gls(rate ~ (pressure + I(pressure^2) + I(pressure^3)) * QB, Dialyzer, 
    method = "ML") 

varying <- list(correlation = alist(AR1_0.771 = corAR1(0.771, form = ~1 | Subject), 
    AR1 = corAR1(), NULL), weights = alist(vp.press = varPower(form = ~pressure), 
    NULL)) 


dd <- dredge(fm1Dial.gls, m.lim = c(1, 2), fixed = ~pressure, varying = varying, 
    extra = "R^2") 

models <- get.models(dd, subset = 1:4)

predict(fm1Dial.gls, se.fit = TRUE, newdata = Dialyzer[1:5, ])

subset(dd, correlation == "AR1_0.771", recalc.delta = TRUE)

ma <- model.avg(models, revised = TRUE)
ms <- model.sel(models)
print(ms, abbr = FALSE)
print(ms, abbr = TRUE)
summary(ma)
predict(ma)[1:10] 


# testing predict replacement:
fm1 <- lme(rate ~ (pressure + I(pressure^2) + I(pressure^3)) * QB, ~1 | Subject, 
    data = Dialyzer)
predict(fm1, newdata = Dialyzer[1:5, ], level = 0, se.fit = TRUE)
 

detach(package:nlme); rm(list=ls())

# TEST glmmML --------------------------------------------------------------------
if (.checkPkg("glmmML")) {
    library("glmmML")
    
    set.seed(100)
    dat <- data.frame(y = rbinom(100, prob = rep(runif(20), rep(5, 20)), size = 1), 
        x = rnorm(100), x2 = rnorm(100), id = factor(rep(1:20, rep(5, 20))))
    
    fm1 <- glmmML(y ~ x * x2, data = dat, cluster = id, x = TRUE)
    dd <- dredge(fm1)
    # mod <- get.models(dd, subset = delta <= 4)
    summary(ma <- model.avg(dd, subset = delta <= 4))
    # vcov(ma)
    coefTable(ma)
    
    detach("package:glmmML")
    rm(list = ls())
} 


# TEST lm ---------------------------------------------------------------------------------
if (.checkPkg("nlme")) {
    data(Orthodont, package = "nlme")
    
    fm1 <- lm(distance ~ Sex * age + age * Sex, data = Orthodont)
    
    dispersion <- function(object) {
        wts <- weights(object)
        if (is.null(wts)) 
            wts <- 1
        sum((wts * resid(object, type = "working")^2)[wts > 0])/df.residual(object)
    }
    
    dd <- dredge(fm1, extra = alist(dispersion))
    gm <- get.models(dd, subset = 1:4)
    ma <- model.avg(gm, revised = F)
    
    vcov(ma)
    summary(ma)
    confint(ma)
    
    predict(ma)
    predict(ma, se.fit = TRUE)
    predict(ma, data.frame(Sex = "Male", age = 8:12))
    
    rm(list = ls())
} 


# TEST glm --------------------------------------------------------------------------------
data(Cement, package = "MuMIn")

nseq <- function(x, len = length(x)) seq(min(x, na.rm = TRUE), 
    max(x, na.rm = T), length = len)
	fm1 = glm(y ~ (X1 + X2 + X3)^2, data = Cement)
	dd <- dredge(fm1)

gm = get.models(dd, subset = 1L:10L)

summary(ma <- model.avg(gm))

vcov(ma)

summary((ma1 = model.avg(dd[1L:10L])))
summary(ma2 <- model.avg(model.sel(dd[1L:10L], rank = "AICc")))
all.equal(ma$avg.model, ma1$avg.mode)

predict(ma) == predict(ma, Cement)
predict(ma, se.fit = T)
predict(ma, lapply(Cement, nseq))


rm(list=ls())
# TEST rlm --------------------------------------------------------------------------------

if (.checkPkg("MASS")) {
library(MASS)
data(Cement, package = "MuMIn")

nseq <- function(x, len=length(x)) seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE),
	length=len)

fm1 <- rlm(y ~X1+X2*X3+X4, data = Cement)
dd <- dredge(fm1, trace=T)
gm <- get.models(dd, subset = 1:10)
ma <- model.avg(gm)
stopifnot(all(predict(ma) == predict(ma, Cement)))
predict(ma, lapply(Cement, nseq, len=30), se.fit=TRUE)
vcov(ma)

rm(list=ls()); #detach(package:MASS)

# TEST multinom --------------------------------------------------------------------------------
if (.checkPkg("nnet")) {
library(nnet); library(MASS)

# Trimmed-down model from example(birthwt)
data(birthwt)

bwt <- with(birthwt, data.frame(
		low = low,
		race = factor(race, labels = c("white", "black", "other")),
		ptd = factor(ptl > 0),
		smoke = (smoke > 0)
		))

options(contrasts = c("contr.treatment", "contr.poly"))
bwt.mu <- multinom(low ~ ., data = bwt)
dd <- dredge(bwt.mu, trace=T)

summary(model.avg(dd[1:5]))
gm <- get.models(dd, subset = 1:5)
ma <- model.avg(gm)

summary(ma)

# predict(ma) // Cannot average factors!

rm(list=ls()); detach(package:nnet)
}}


# TEST gam --------------------------------------------------------------------------------
if (.checkPkg("mgcv")) {

suppressPackageStartupMessages(library(mgcv))
RNGkind("Mersenne")
set.seed(0) ## simulate some data...
dat <- gamSim(1, n = 400, dist = "binary", scale = 2)
#gam1 <- gam(y~s(x0)+s(x1)+s(x2)+s(x3), data=dat)

ops <- options(warn = -1)

gam1 <- gam(y ~ s(x0) + s(x1) + s(x2) +  s(x3) + (x1+x2+x3)^2,
	data = dat, method = "GCV.Cp", family = binomial)



dd <- dredge(gam1, subset=!`s(x0)` & (!`s(x1)` | !x1) & (!`s(x2)` |
	!x2) & (!`s(x3)` | !x3), fixed = "x1")
	
gm <- get.models(dd, cumsum(weight) <= .95)
ma <- model.avg(gm)

summary(ma)

predict(ma, dat[1:10, ], se.fit=T, type = "link")

predict(ma, dat[1:10, ], se.fit=T, type = "response")
predict(ma, dat[1:10, ], se.fit=T, type = "link", backtransform = TRUE)

options(ops)

rm(list=ls()); detach(package:mgcv)
}

# TEST spautolm ---------------------------------------------------------------------------

# if (require("foreign") && require("spdep"))
if (.checkPkg(c("foreign", "spdep")))
if(!is.null(tryCatch(suppressPackageStartupMessages(library(spdep)), error = function(e) NULL))) {

suppressMessages(example(NY_data, echo = FALSE))

# method argument changed in spdep 0.5.53 from "full" to "eigen"
method1 <- if(packageVersion("spdep") < '0.5.53') "full" else "eigen"


fm1.spautolm <- spautolm(Z ~ PEXPOSURE * PCTAGE65P + PCTOWNHOME,
 data = nydata, listw = listw_NY, family = "SAR", method = method1, verbose = FALSE)

options(warn=1)
dd <- dredge(fm1.spautolm, m.lim=c(0,1), fixed = ~PEXPOSURE,
	varying = list(
		family = list("CAR", "SAR"),
		method=list("Matrix_J", method1)
	), trace=FALSE)
options(warn=0)


#dd <- dredge(fm1.spautolm, m.lim=c(0,3), fixed=~PEXPOSURE)
gm <- get.models(dd, cumsum(weight) <= .99)
ma <- model.avg(gm)
summary(ma)
# signif(resid(ma), 5)[1:10]

rm(list=ls())

# TEST spautolm ---------------------------------------------------------------------------
#suppressPackageStartupMessages(library(spdep))
data(oldcol)

fm1.sarlm <- errorsarlm(CRIME ~ INC * HOVAL * OPEN, data = COL.OLD,
 listw = nb2listw(COL.nb, style = "W"), method = "eigen", quiet = TRUE)

 
dd <- dredge(fm1.sarlm)

gm <- get.models(dd, cumsum(weight) <= .98)
ma <- model.avg(gm)


#fm <- fm1.sarlm
#avgpred(ma)
#avgpred(ma, newdata = COL.OLD, listw = nb2listw(COL.nb, style = "W"))
#avgpred(ma, newdata = COL.OLD[1:10, ])
#
#std_predict(fm, newdata = COL.OLD, listw = nb2listw(COL.nb, style = "W"))
#predict(fm, newdata = COL.OLD, listw = nb2listw(COL.nb, style = "W"))

stopifnot(isTRUE(all.equal(coefTable(ma), coefTable(model.avg(dd, cumsum(weight) <= .98)))))

summary(ma)

predict(ma)[1:10]

rm(list=ls()); detach(package:spdep)

} # library(spdep)

# TEST glm.nb ---------------------------------------------------------------------------
if (.checkPkg("MASS")) {

require("MuMIn")
options(na.action = na.fail)
require("MASS")

quine.nb1 <- glm.nb(Days ~ 0 + Sex/(Age + Eth*Lrn), data = quine)
#quine.nb1 <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine)

ms <- dredge(quine.nb1)

models <- get.models(ms, subset = TRUE)
models <- get.models(ms, subset = NA)

summary(model.avg(models))



#dredge(quine.nb1) # OK
#dredge(quine.nb1x = NA) # OK
#dredge(quine.nb1) # OK
dredge(quine.nb1) # OK
#dredge(quine.nb1) # Right, should be the same as above
ma <- model.avg(dredge(quine.nb1), subset = cumsum(weight)<=.9999)

# Cannot predict with this 'averaging'
#pred <- predict(ma, se=TRUE)

#pred <- cbind(pred$fit, pred$fit - (2 * pred$se.fit), pred$fit + (2 * pred$se.fit))
#matplot(pred, type="l")
#matplot(family(quine.nb1)$linkinv(pred), type="l")

rm(list=ls()); #detach(package:MASS)
}

# TEST quasibinomial -----------------------------------------------------------

budworm <- data.frame(ldose = rep(0:5, 2), numdead = c(1, 4, 9, 13, 18, 20, 0,
	2, 6, 10, 12, 16), sex = factor(rep(c("M", "F"), c(6, 6))))
budworm$SF = cbind(numdead = budworm$numdead, numalive = 20 - budworm$numdead)

qbinomial <- function(...) {
	res <- quasibinomial(...)
	res$aic <- binomial(...)$aic
	res
}
budworm.qqlg <- glm(SF ~ sex*ldose + sex*I(ldose^2), data = budworm, family = qbinomial)
#budworm.qlg <- glm(SF ~ sex*ldose + sex*I(ldose^2), data = budworm, family = quasibinomial)
budworm.lg <- glm(SF ~ sex*ldose + sex*I(ldose^2), data = budworm, family = binomial)

#R2(budworm.lg)
r.squaredLR(budworm.lg)
r.squaredGLMM(budworm.lg)

dd <- dredge(budworm.lg, rank = "QAIC", chat = summary(budworm.lg)$dispersion)
#dd <- dredge(budworm.lg) # should be the same
mod <- get.models(dd, subset = NA)

# Note: this works:
# ma <- model.avg(mod)
# but this will not: ('rank' attribute passed from 'dredge' is lost)
# ma <- model.avg(mod)
# so, need to supply them
ma <- model.avg(mod[1:5], rank = "QAICc", rank.args = list(chat = 0.403111))

rm(list=ls())

# TEST polr {MASS} -------------------------------------------------------------
#if (.checkPkg("MASS")) {
#library(MASS)
#library(MuMIn)
#options(contrasts = c("contr.treatment", "contr.poly"))
#house.plr <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)

#dd <- dredge(house.plr)
#mod <- get.models(dd, 1:6)
#model.avg(mod)
#}


# TEST coxph -------------------------------------------------------------------

library(survival)

bladder1 <- bladder[bladder$enum < 5, ]

fmcph <- coxph(Surv(stop, event) ~ (rx + size + number) * strata(enum) +  cluster(id), bladder1)

r.squared.coxph <- function(object, ...) {
	logtest <- -2 * (object$loglik[1L] - object$loglik[2L])
	c(rsq = 1 - exp(-logtest/object$n), maxrsq = 1 - exp(2 * object$loglik[1L]/object$n))
}

ms <- dredge(fmcph, fixed=c("cluster(id)", "strata(enum)"), extra = list(R2="r.squared.coxph"))

fits <- get.models(ms, delta < 5)
summary(model.avg(fits))

fmsrvrg <- survreg(Surv(futime, fustat) ~ ecog.ps + rx, ovarian, dist='weibull',
    scale = 1)

summary(model.avg(dredge(fmsrvrg), delta  < 4))

fmsrvrg2 <- survreg(Surv(futime, fustat) ~ ecog.ps + rx, ovarian, dist='weibull')

fmsrvrg3 <- survreg(Surv(time, status) ~ ph.ecog + age + strata(sex), lung,
	  na.action = "na.omit")


coefTable(fmsrvrg)
coefTable(fmsrvrg2)
coefTable(fmsrvrg3)


rm(list=ls())
detach(package:survival)

# END TESTS
