### R code from vignette source 'rlmer.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: initialize
###################################################

#####################################################
#####################################################
## set this to true to update fitted objects cache ##
UPDATE_CACHE <- FALSE
#####################################################
#####################################################

options(width=80,
        str=strOptions(strict.width = "wrap", vec.len=3),
        xtable.size="\\small")
if (!file.exists("figs")) dir.create("figs")
rlmerDoc <- system.file("doc", package="robustlmm")
require(ggplot2)
require(reshape2)
source(file.path(rlmerDoc, "ggplot.theme.R"))
require(xtable)
require(robustlmm)
require(robustbase)
require(lme4)
source(file.path(rlmerDoc, "plots.R"))
st <- system.time
if(packageVersion("robustbase") >= "0.9-8")  {
    lqqPsi <- psiFuncCached(rho = function(x, cc) Mpsi(x, cc, "lqq", -1),
                            psi = function(x, cc) Mpsi(x, cc, "lqq", 0),
                            Dpsi = function(x, cc) Mpsi(x, cc, "lqq", 1),
                            wgt = function(x, cc) Mwgt(x, cc, "lqq"),
                            Dwgt = function(x, cc) 
                            (Mpsi(x, cc, "lqq", 1) - Mwgt(x, cc, "lqq"))/x,
                            name = "lqq",
                            cc = c(-0.5, 1.5, 0.95, NA))
    
    bisquarePsi <- psiFuncCached(rho = function(x, k) Mpsi(x, k, "biweight", -1),
                                 psi = function(x, k) Mpsi(x, k, "biweight", 0),
                                 Dpsi = function(x, k) Mpsi(x, k, "biweight", 1),
                                 wgt = function(x, k) (1 - (x/k)^2)^2*(abs(x) <= k),
                                 Dwgt = function(x, k) (-(4*(1-(x/k)^2))*x/k^2)*(abs(x) <= k),
                                 name = "bisquare",
                                 k = 4.68)
} else {
    lqqPsi <- psiFuncCached(rho = function(x, cc) robustbase:::lmrob.psifun(x, cc, "lqq", -1),
                            psi = function(x, cc) robustbase:::lmrob.psifun(x, cc, "lqq", 0),
                            Dpsi = function(x, cc) robustbase:::lmrob.psifun(x, cc, "lqq", 1),
                            wgt = function(x, cc) robustbase:::lmrob.wgtfun(x, cc, "lqq"),
                            Dwgt = function(x, cc) (robustbase:::lmrob.psifun(x, cc, "lqq", 1) -
                                robustbase:::lmrob.wgtfun(x, cc, "lqq"))/x,
                            name = "lqq",
                            cc = c(-0.5, 1.5, 0.95, NA))
    
    bisquarePsi <- psiFuncCached(rho = function(x, k) tukeyPsi1(x, k, -1),
                                 psi = function(x, k) tukeyPsi1(x, k, 0),
                                 Dpsi = function(x, k) tukeyPsi1(x, k, 1),
                                 wgt = function(x, k) (1 - (x/k)^2)^2*(abs(x) <= k),
                                 Dwgt = function(x, k) (-(4*(1-(x/k)^2))*x/k^2)*(abs(x) <= k),
                                 name = "bisquare",
                                 k = 4.68)
}
## make the functions wgt.e and wgt.b accessible again...
wgt.e <- robustlmm:::wgt.e
wgt.b <- robustlmm:::wgt.b


###################################################
### code chunk number 2: penicillin-setup
###################################################
data(Penicillin)
## setup datasets
Penicillin <- within(Penicillin, {
    plate <- reorder(plate, diameter)
    attr(plate, "scores") <- NULL
})


###################################################
### code chunk number 3: penicillin-raw
###################################################
print(ggplot(Penicillin, aes(plate, diameter, color = sample)) +
      geom_point() + geom_line(aes(as.numeric(plate))) +
      scale_colour_brewer("Sample", palette="Dark2") +
      scale_y_continuous(breaks=c(18,20,22,24,26)) +
      xlab("Plate") + ylab("Diameter growth inhibition zone (mm)") +
      theme(legend.position = "bottom", legend.box = "horizontal"))


###################################################
### code chunk number 4: penicillin-str
###################################################
str(Penicillin)


###################################################
### code chunk number 5: penicillin-lmer
###################################################
st(classical <- lmer(diameter ~ 1 + (1|plate) + (1|sample),
                     Penicillin))


###################################################
### code chunk number 6: penicillin-rlmer
###################################################
st(robust <- rlmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin,
                   rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
                   rho.sigma.b = psi2propII(smoothPsi, k = 2.28)))
summary(robust)


###################################################
### code chunk number 7: penicillin2
###################################################
st(robust2 <- rlmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin,
                    rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
                    rho.b = list(smoothPsi, cPsi),
                    rho.sigma.b = list(psi2propII(smoothPsi, k = 2.28),
                                       cPsi)))


###################################################
### code chunk number 8: penicillin-cmp
###################################################
print(xtable(compare(classical, robust, robust2, show.rho.functions=FALSE),
             caption="Comparison table of the fitted models for the Penicillin example.",
             label="tab:cmpPenicillin"))


###################################################
### code chunk number 9: penicillin-ta
###################################################
plots <- plot(robust)
lower <- floor(min(getME(robust, "w_e"), getME(robust, "w_b_vector"))*100)/100
plots[[1]] + scale_colour_gradient(limits=c(lower,1)) +
    theme(legend.position = "none")


###################################################
### code chunk number 10: penicillin-qq-resid
###################################################
plots[[2]] + scale_colour_gradient(limits=c(lower,1)) +
    theme(legend.position = "none")


###################################################
### code chunk number 11: penicillin-qq-ranef
###################################################
plots[[3]] +
    scale_colour_gradient("robustness weights", limits=c(lower,1)) +
    theme(legend.position = "bottom", legend.box = "horizontal")


###################################################
### code chunk number 12: penicillin-robustness-weights
###################################################
tmp <- cbind(Penicillin, wgt.e = getME(robust, "w_e"))
print(ggplot(tmp, aes(plate, diameter, color = sample)) +
      geom_point(aes(size=1/wgt.e)) + geom_line(aes(as.numeric(plate))) +
      scale_colour_brewer("Sample", palette="Dark2") +
      scale_size_continuous(expression(w[e]),breaks=c(1,1.5,2,2.5),labels=c(1,0.66,0.5,0.33)) +
      scale_y_continuous(breaks=c(18,20,22,24,26)) +
      xlab("Plate") + ylab("Diameter growth inhibition zone (mm)") +
      theme(legend.position = "bottom", legend.box = "horizontal"))


###################################################
### code chunk number 13: sleepstudy-setup
###################################################
source(file.path(rlmerDoc, "sleepstudy.R"))


###################################################
### code chunk number 14: sleepstudy-raw
###################################################
print(ggplot(sleepstudy, aes(Days, Reaction)) +
      stat_smooth(method=function(formula, ..., weights)
                  lmrob(formula, ..., setting="KS2011"),
                  se=FALSE) +
      ## stat_smooth(method="loess", span=1.5, se=FALSE) +
      geom_point() +
      scale_x_continuous(breaks=c(0,2,4,6,8)) +
      facet_wrap(~ Subject, nrow=3) +
      xlab("Days of sleep deprivation") +
      ylab("Average reaction time (ms)"))


###################################################
### code chunk number 15: sleepstudy-str
###################################################
str(sleepstudy)


###################################################
### code chunk number 16: sleepstudy (eval = FALSE)
###################################################
## st(classical <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))


###################################################
### code chunk number 17: sleepstudy3 (eval = FALSE)
###################################################
## st(robust <-
##    rlmer(Reaction ~ Days + (Days|Subject), sleepstudy,
##          rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
##          rho.sigma.b = chgDefaults(smoothPsi, k = 5.11, s=10)))
## summary(robust)


###################################################
### code chunk number 18: smoothedHuber
###################################################
xs <- seq.int(0, 3, length.out=100)
data <- data.frame(x = xs,
                   Huber = huberPsi@psi(xs),
                   `Smoothed` = smoothPsi@psi(xs))
print(ggplot(melt(data, 1),
             aes(x, value, color=variable)) + geom_line() +
      scale_colour_hue(expression(paste(psi, "-function"))) +
      ylab(expression(psi(x))) +
      theme(legend.position = "bottom", legend.box = "horizontal"))



###################################################
### code chunk number 19: efficiency-table (eval = FALSE)
###################################################
## require(robustlmm)
## 
## ## location
## avarBeta <- function(psi) psi@Epsi2()/psi@EDpsi()^2
## effBeta <- function(psi) 1 / avarBeta(psi)
## findEffBeta <- function(psi, eff) {
##     tDefs <- psi@tDefs
##     if (length(tDefs) == 0) return(NULL)
##     lavar <- function(c) {
##         tDefs[1] <- c
##         lpsi <- do.call(chgDefaults, c(list(psi), tDefs))
##         effBeta(lpsi) - eff
##     }
##     ret <- try(uniroot(lavar, tDefs[1] * c(0.25, 5))$root)
##     ## fail gracefully
##     if (is(ret, "try-error")) NA else ret
## }
## 
## ## scale
## E <- function(expr, ...) {
##     sexpr <- substitute(expr)
##     args <- list(...)
##     if (is.name(sexpr)) {
##         fcall <- paste(sexpr, "(x)")
##         expr <- parse(text = fcall)
##     }
##     else {
##         if (!(is.call(sexpr) && match("x", all.vars(sexpr), nomatch = 0L)))
##             stop("'expr' must be a function or an expression containing 'x'")
##         expr <- sexpr
##     }
##     ifun <- function(x)
##         eval(expr, envir = c(list(x=x), args), enclos=parent.frame())*dnorm(x)
##     integrate(ifun, -Inf, Inf)$value
## }
## 
## ## kappa = E[w(x)x^2] / E[w(x)]
## kappa <- function(k) E(w(x, k)*x^2, k=k) / E(w(x, k), k=k)
## ## Psi(x,sigma,k) = w(x/sigma, k)*((x/sigma)^2 - kappa(k))
## Psi <- function(x, sigma, k) w(x/sigma,k)*((x/sigma)^2 - kappa(k))
## ## Dpsi(x,sigma,k) = d/dsigma Psi(x,sigma,k)
## DPsi <- function(x, sigma, k) {
##     xs <- x/sigma
##     -1/sigma*(Dw(xs,k)*((xs)^3 - kappa(k)*(xs)) + 2*w(xs,k)*xs^2)
## }
## ## A = E[Psi(x,1,k)^2]
## A <- function(k) E(Psi(x,1,k)^2,k=k)
## ## B = E[d/dsigma Psi
## B <- function(k) E(DPsi(x,1,k), k=k)
## ## effSigma = B / A^2 / 2
## effSigma <- Vectorize(function(k) B(k)^2 / A(k) / 2)
## ## find eff:
## findEffSigma <- function(eff, lower=0.4, upper=6) {
##     lfun <- function(k) effSigma(k) - eff
##     uniroot(lfun, lower=lower, upper=upper)$root
## }
## 
## ## compute tuning constants to give the desired efficiency
## require(xtable)
## effs <- c(0.80, 0.85, 0.90, 0.95)
## tbl <- matrix("", nrow=length(effs), ncol=4)
## tbl[,1] <- sprintf("%.2f", effs)
## colnames(tbl) <- c("efficiency",
##                    "$k$ for $\\hat\\mu$",
##                    "$k$ for $\\hat\\sigma_\\text{\\eqref{eq:sigmaDAS}}$",
##                    "$k$ for $\\hat\\sigma_\\text{\\eqref{eq:sigmaDAS}}$, Prop. II")
## ## location and wls case
## ## w and Dw for wls case
## w <- huberPsi@wgt
## Dw <- huberPsi@Dwgt
## for (i in seq_along(effs)) {
##     leff <- effs[i]
##     tbl[i,2] <- round(findEffBeta(huberPsi, leff), 2)
##     tbl[i,3] <- round(findEffSigma(leff), 2)
## }
## ## Prop II case
## huberProp2 <- psi2propII(huberPsi)
## w <- huberProp2@wgt
## Dw <- huberProp2@Dwgt
## for (i in seq_along(effs)) {
##     leff <- effs[i]
##     tbl[i,4] <- round(findEffSigma(leff), 2)
## }
## 
## ## print results
## print(xtable(tbl, align="lcccc"),
##       file="tuning-constants-table.tex",include.rownames=FALSE,
##       floating=FALSE,sanitize.text.function=function(x){x})
## 


###################################################
### code chunk number 20: efficiency-tables (eval = FALSE)
###################################################
## ## dist functions
## dist <- function(b, kappa=kappa) {
##     if (!is.matrix(b)) return(b*b) ## assume s==1
##     s <- ncol(b)
##     if (s == 1) return(drop(b*b))
##     ## else: just square and sum
##     rowSums(b*b) - kappa*s
## }
## dist2 <- function(sb2, s=1, kappa=kappa) {
##     if (s==1) return(sb2)
##     return(sb2 - kappa*s)
## }
## kappa.nondiag <- function(wgt, s) {
##     if (s == 1) {
##         tfun <- function(v) wgt(dist2(v))*v*dchisq(v,s)
##         tfun2 <- function(v) wgt(dist2(v))*dchisq(v,s)
##         integrate(tfun, 0, Inf)$value / integrate(tfun2, 0, Inf)$value
##     } else {
##         tfun <- function(v, kappa) (v-s*kappa)*wgt(v - s*kappa)*dchisq(v, s)
##         tfun2 <- function(kappa) integrate(tfun, 0, Inf, kappa=kappa)$value
##         uniroot(tfun2, c(0, 1))$root
##     }
## }
## ## asymptotic efficiency functions
## dEta <- Vectorize(function(wgt, s) {
##     ## kappa <- kappa.nondiag(wgt, s)
##     DEtaQ <- 1/(1 + 2/s) * integrate(function(v) (v/s)^2*wgt(v)^2*dchisq(v,s), 0, Inf)$value
##     DEtaM <- 1/(1 + 2/s) * integrate(function(v) (v/s)^2*wgt(v)*dchisq(v,s), 0, Inf)$value
##     DEtaM^2 / DEtaQ
## }, "s")
## dTau <- Vectorize(function(wgt, s) {
##     kappa <- kappa.nondiag(wgt,s)
##     uTau <- function(v) dist2(v,s,kappa)*wgt(dist2(v,s,kappa))
##     DTauQ <- 1/(2*s) * integrate(function(v) uTau(v)^2*dchisq(v,s), 0, Inf)$value
##     DTauM <- 1/(2*s) * integrate(function(v) uTau(v)*(v-s)*dchisq(v,s), 0, Inf)$value
##     DTauM^2 / DTauQ
## }, "s")
## dMu <- Vectorize(function(wgt, s) {
##     DMuQ <- 1/s * integrate(function(v) v*wgt(sqrt(v))^2*dchisq(v,s), 0, Inf)$value
##     DMuM <- 1/s * integrate(function(v) v*wgt(sqrt(v))*dchisq(v,s), 0, Inf)$value
##     DMuM^2 / DMuQ
## }, "s")
## ## function to find corresponding efficiency
## findCorrespEff.nondiag <- function(rho, s, effFun, eff=effBeta(rho), lower=0.1) {
##     cat("Eff for beta:", eff, "\n")
##     lfun <- function(k) {
##         lrho <- chgDefaults(rho, k = k)
##         effFun(lrho@wgt, s) - eff
##     }
##     uniroot(lfun, c(lower, 15))$root
## }
## ## need a different function for lqq
## findCorrespEff.lqq <- function(s, effFun, eff, lower=1) {
##     if (s == 1) {
##         lfun <- function(k) {
##             lrho <- chgDefaults(lqqPsi, cc = c(1.5*k, k, 1.5))
##             asymptEff1(lrho@wgt, lrho@Dwgt) - eff
##         }
##     } else {
##         lfun <- function(k) {
##             lrho <- chgDefaults(lqqPsi, cc = c(1.5*k, k, 1.5))
##             effFun(lrho@wgt, s) - eff
##         }
##     }
##     k. <- uniroot(lfun, c(lower, 15))$root
##     c(1.5*k., k., 1.5)
## }
## 
## ## compute table for huber functions...
## ss <- c(2:7)
## tbl <- matrix("NA", nrow=4, ncol=length(ss))
## tbl[1,] <- as.character(ss)
## for (i in seq_along(ss)) {
##     s <- ss[i]
##     tbl[2,i] <- format(findCorrespEff.nondiag(huberPsi, s, dEta, 0.95), digits=3)
##     tbl[3,i] <- format(findCorrespEff.nondiag(huberPsi, s, dTau, 0.95), digits=3)
##     tbl[4,i] <- format(findCorrespEff.nondiag(huberPsi, s, dMu, 0.95), digits=3)
## }
## rownames(tbl) <- c("$s$", "$b_\\eta$", "$b_\\tau$", "$b_\\mu$")
## print(xtable(tbl,
##              caption="Tuning parameters for the optimal $B$-estimator to yield $95\\%$ efficiency, non-diagonal case. For the Huber $\\psi$-function.",
##              label="tab:effBOptimal"),
##       hline.after=c(0,1,4),
##       include.colnames=FALSE,
##       sanitize.text.function=identity,
##       file="efficiency-table-B-optimal.tex")
## 
## ## compute loc-scale table for lqq functions...
## effs <- c(0.80, 0.85, 0.90, 0.95)
## tbl <- matrix("", nrow=length(effs), ncol=3)
## tbl[,1] <- sprintf("%.2f", effs)
## f <- function(cc) sprintf("($%s$)", paste(format(cc[-3], digits=3),collapse="$,$"))
## for (i in seq_along(effs)) {
##     leff <- effs[i]
##     tbl[i,2] <- f(robustbase:::lmrob.lqq.findc(c(-0.5, 1.5, leff, NA)))
##     tbl[i,3] <- f(findCorrespEff.lqq(1, eff=leff, lower=0.5))
## }
## colnames(tbl) <- c("efficiency", "$cc$ for $\\hat\\mu$",
##                    "$cc$ for $\\hat\\sigma_\\text{\\eqref{eq:sigmaDAS}}$")
## print(xtable(tbl,
##              caption="Tuning parameters for lqq $\\psi$-function for the location and scale estimates such that they reach the given asymptotic efficiency. The third parameter is always taken to be $1.5$.",
##              label="tab:eqqLqqLocationScale"),
##       include.rownames=FALSE,
##       sanitize.text.function=identity,
##       file="efficiency-table-lqq-locationScale.tex")
## 
## ## compute non-diag table for lqq functions...
## ss <- c(2:6)
## tbl <- matrix("", nrow=4, ncol=length(ss))
## tbl[1,] <- as.character(ss)
## for (i in seq_along(ss)) {
##     s <- ss[i]
##     tbl[2,i] <- f(findCorrespEff.lqq(s, dEta, 0.95))
##     tbl[3,i] <- f(findCorrespEff.lqq(s, dTau, 0.95))
##     tbl[4,i] <- f(findCorrespEff.lqq(s, dMu, 0.95))
## }
## rownames(tbl) <- c("$s$", "$cc_\\eta$", "$cc_\\tau$", "$cc_\\mu$")
## print(xtable(tbl, digits=2,
##              caption="Tuning parameters for the lqq weight function to yield $95\\%$ efficiency, non-diagonal case. The third parameter is always taken to be $1.5$.",
##              label="tab:effLqq"),
##       hline.after=c(0,1,4),
##       include.colnames=FALSE,
##       sanitize.text.function=identity,
##       file="efficiency-table-lqq.tex")


