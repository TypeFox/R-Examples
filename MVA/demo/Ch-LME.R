### R code from vignette source 'Ch-LME.Rnw'

###################################################
### code chunk number 1: setup
###################################################
library("MVA")
set.seed(280875)
library("lattice")
lattice.options(default.theme =
    function()
        standard.theme("pdf", color = FALSE))

if (file.exists("deparse.R")) {
    if (!file.exists("figs")) dir.create("figs")
    source("deparse.R")
    options(prompt = "R> ", continue = "+  ", width = 64,
        digits = 4, show.signif.stars = FALSE, useFancyQuotes = FALSE)

    options(SweaveHooks = list(onefig =   function() {par(mfrow = c(1,1))},
                           twofig =   function() {par(mfrow = c(1,2))},                           
                           figtwo =   function() {par(mfrow = c(2,1))},                           
                           threefig = function() {par(mfrow = c(1,3))},
                           figthree = function() {par(mfrow = c(3,1))},
                           fourfig =  function() {par(mfrow = c(2,2))},
                           sixfig =   function() {par(mfrow = c(3,2))},
                           nomar = function() par("mai" = c(0, 0, 0, 0))))
}


###################################################
### code chunk number 2: setup
###################################################
library("nlme")


###################################################
### code chunk number 3: ch:LME:ex
###################################################
exd <- data.frame(ID = factor(1:6), Group = gl(2, 3), 
                  matrix(rpois(6 * 4, lambda = 10), nrow = 6))
colnames(exd)[-(1:2)] <- paste("Day", c(1, 2, 5, 7), sep = ".")
ex_wide <- exd


###################################################
### code chunk number 4: ch:LME:wide
###################################################
ex_wide


###################################################
### code chunk number 5: ch:LME:wide
###################################################
reshape(ex_wide, direction = "long", idvar = "ID", 
        varying = colnames(ex_wide)[-(1:2)])


###################################################
### code chunk number 6: ch:LME:timber:tab
###################################################
"timber" <-
matrix(c(0., 0., 0., 0., 0., 0., 0., 0., 2.3799999999999999, 2.6899999999999999, 2.8500000000000001, 2.46,
        2.9700000000000002, 3.96, 3.1699999999999999, 3.3599999999999999, 4.3399999999999999, 4.75,
        4.8899999999999997, 4.2800000000000002, 4.6799999999999997, 6.46, 5.3300000000000001, 5.4500000000000002,
        6.6399999999999997, 7.04, 6.6100000000000003, 5.8799999999999999, 6.6600000000000001, 8.1400000000000006,
        7.1399999999999997, 7.0800000000000001, 8.0500000000000007, 9.1999999999999993, 8.0899999999999999,
        7.4299999999999997, 8.1099999999999994, 9.3499999999999996, 8.2899999999999991, 8.3200000000000003,
        9.7799999999999994, 10.94, 9.7200000000000006, 8.3200000000000003, 9.6400000000000006, 10.720000000000001,
        9.8599999999999994, 9.9100000000000001, 10.970000000000001, 12.23, 11.029999999999999, 9.9199999999999999,
        11.06, 11.84, 11.07, 11.06, 12.050000000000001, 13.19, 12.140000000000001, 11.1, 12.25, 12.85,
        12.130000000000001, 12.210000000000001, 12.98, 14.08, 13.18, 12.23, 13.35, 13.83, 13.15, 13.16, 13.94,
        14.66, 14.119999999999999, 13.24, 14.539999999999999, 14.85, 14.09, 14.050000000000001, 14.74,
        15.369999999999999, 15.09, 14.19, 15.529999999999999, 15.789999999999999, 15.109999999999999,
        14.960000000000001, 16.129999999999999, 16.890000000000001, 16.68, 16.07, 17.379999999999999,
        17.390000000000001, 16.690000000000001, 16.239999999999998, 17.98, 17.780000000000001, 17.940000000000001,
        17.43, 18.760000000000002, 18.440000000000001, 17.690000000000001, 17.34, 19.52, 18.41, 18.219999999999999,
        18.359999999999999, 19.809999999999999, 19.460000000000001, 18.710000000000001, 18.23, 19.969999999999999,
        18.969999999999999, 19.399999999999999, 18.93, 20.620000000000001, 20.050000000000001, 19.539999999999999,
        18.870000000000001)
, nrow = 8, ncol = 15)

slippage <- c((0:10)/10, seq(from = 1.2, to = 1.8, by = 0.2))

colnames(timber) <- paste("s", slippage, sep = "")
timber <- as.data.frame(timber)
timber$specimen <- factor(paste("spec", 1:nrow(timber), sep = ""))

timber.dat <- reshape(timber, direction = "long", idvar = "specimen",
        varying = matrix(colnames(timber)[1:15], nr = 1),
        timevar = "slippage")
names(timber.dat)[3] <- "loads"
timber.dat$slippage <- slippage[timber.dat$slippage]
timber <- timber.dat

toLatex(HSAURtable(timber), pcol = 2,
    caption = paste("Data giving loads needed for a given slippage",
                    "in eight specimens of timber, with data in ``long'' form."),
    label = "ch:LME:timber:tab", rownames = FALSE)


###################################################
### code chunk number 7: ch:LME:plasma:tab
###################################################
"plasma" <- 
matrix(c(4.2999999999999998, 3.7000000000000002, 4., 3.6000000000000001, 4.0999999999999996, 3.7999999999999998,
        3.7999999999999998, 4.4000000000000004, 5., 3.7000000000000002, 3.7000000000000002, 4.4000000000000004,
        4.7000000000000002, 4.2999999999999998, 5., 4.5999999999999996, 4.2999999999999998, 3.1000000000000001,
        4.7999999999999998, 3.7000000000000002, 5.4000000000000004, 3., 4.9000000000000004, 4.7999999999999998,
        4.4000000000000004, 4.9000000000000004, 5.0999999999999996, 4.7999999999999998, 4.2000000000000002,
        6.5999999999999996, 3.6000000000000001, 4.5, 4.5999999999999996, 3.2999999999999998, 2.6000000000000001,
        4.0999999999999996, 3., 3.7999999999999998, 2.2000000000000002, 3., 3.8999999999999999, 4.,
        3.1000000000000001, 2.6000000000000001, 3.7000000000000002, 3.1000000000000001, 3.2999999999999998,
        4.9000000000000004, 4.4000000000000004, 3.8999999999999999, 3.1000000000000001, 5., 3.1000000000000001,
        4.7000000000000002, 2.5, 5., 4.2999999999999998, 4.2000000000000002, 4.2999999999999998, 4.0999999999999996,
        4.5999999999999996, 3.5, 6.0999999999999996, 3.3999999999999999, 4., 4.4000000000000004, 3.,
        2.6000000000000001, 3.1000000000000001, 2.2000000000000002, 2.1000000000000001, 2., 2.3999999999999999,
        2.7999999999999998, 3.3999999999999999, 2.8999999999999999, 2.6000000000000001, 3.1000000000000001,
        3.2000000000000002, 3., 4.0999999999999996, 3.8999999999999999, 3.1000000000000001, 3.2999999999999998,
        2.8999999999999999, 3.2999999999999998, 3.8999999999999999, 2.2999999999999998, 4.0999999999999996,
        4.7000000000000002, 4.2000000000000002, 4., 4.5999999999999996, 4.5999999999999996, 3.7999999999999998,
        5.2000000000000002, 3.1000000000000001, 3.7000000000000002, 3.7999999999999998, 2.6000000000000001,
        1.8999999999999999, 2.2999999999999998, 2.7999999999999998, 3., 2.6000000000000001, 2.5, 2.1000000000000001,
        3.3999999999999999, 2.2000000000000002, 2.2999999999999998, 3.2000000000000002, 3.2999999999999998,
        2.6000000000000001, 3.7000000000000002, 3.8999999999999999, 3.1000000000000001, 2.6000000000000001, 
        2.7999999999999998, 2.7999999999999998, 4.0999999999999996, 2.2000000000000002, 3.7000000000000002, 
        4.5999999999999996, 3.3999999999999999, 4., 4.0999999999999996, 4.4000000000000004, 3.6000000000000001, 
        4.0999999999999996, 2.7999999999999998, 3.2999999999999998, 3.7999999999999998, 2.2000000000000002, 
        2.8999999999999999, 2.8999999999999999, 2.8999999999999999, 3.6000000000000001, 3.7999999999999998,
        3.1000000000000001, 3.6000000000000001, 3.2999999999999998, 1.5, 2.8999999999999999, 3.7000000000000002,
        3.2000000000000002, 2.2000000000000002, 3.7000000000000002, 3.7000000000000002, 3.1000000000000001,
        2.6000000000000001, 2.2000000000000002, 2.8999999999999999, 2.7999999999999998, 2.1000000000000001,
        3.7000000000000002, 4.7000000000000002, 3.5, 3.2999999999999998, 3.3999999999999999, 4.0999999999999996,
        3.2999999999999998, 4.2999999999999998, 2.1000000000000001, 2.3999999999999999, 3.7999999999999998, 2.5,
        3.2000000000000002, 3.1000000000000001, 3.8999999999999999, 3.3999999999999999, 3.6000000000000001,
        3.3999999999999999, 3.7999999999999998, 3.6000000000000001, 2.2999999999999998, 2.2000000000000002,
        4.2999999999999998, 4.2000000000000002, 2.5, 4.0999999999999996, 4.2000000000000002, 3.1000000000000001,
        1.8999999999999999, 3.1000000000000001, 3.6000000000000001, 3.7000000000000002, 2.6000000000000001,
        4.0999999999999996, 3.7000000000000002, 3.3999999999999999, 4.0999999999999996, 4.2000000000000002, 4.,
        3.1000000000000001, 3.7999999999999998, 2.3999999999999999, 2.2999999999999998, 3.6000000000000001,
        3.3999999999999999, 3.1000000000000001, 3.8999999999999999, 3.7999999999999998, 3.6000000000000001, 3.,
        3.5, 4., 4., 2.7000000000000002, 3.1000000000000001, 3.8999999999999999, 3.7000000000000002,
        2.3999999999999999, 4.7000000000000002, 4.7999999999999998, 3.6000000000000001, 2.2999999999999998, 3.5,
        4.2999999999999998, 3.5, 3.2000000000000002, 4.7000000000000002, 3.6000000000000001, 3.7999999999999998,
        4.2000000000000002, 4.4000000000000004, 3.7999999999999998, 3.5, 4.2000000000000002, 2.5,
        3.1000000000000001, 3.7999999999999998, 4.4000000000000004, 3.8999999999999999, 4., 4., 3.7000000000000002,
        3.5, 3.7000000000000002, 3.8999999999999999, 4.2999999999999998, 2.7999999999999998, 3.8999999999999999,
        4.7999999999999998, 4.2999999999999998, 3.3999999999999999, 4.9000000000000004, 5., 4., 2.7000000000000002,
        3.6000000000000001, 4.4000000000000004, 3.7000000000000002, 3.5, 4.9000000000000004, 3.8999999999999999,
        4., 4.2999999999999998, 4.9000000000000004, 3.7999999999999998, 3.8999999999999999, 4.7999999999999998, 
        3.5, 3.2999999999999998, 3.7999999999999998)
, nrow = 33, ncol = 8
,  dimnames = list(character(0)
, c("T0.0", "T0.5", "T1.0", "T1.5", "T2.0", "T2.5", "T3.0", "T4.0")
)
)

time <- c(0, 0.5, 1, 1.5, 2, 3, 4)
plasma <- as.data.frame(plasma)
plasma$Subject <- factor(paste("id", 
    formatC(1:nrow(plasma), format = "g", width = 2, flag = "0"), sep = ""))
plasma$group <- factor(c(rep("control", 20), rep("obese", 13)))

plasma <- reshape(plasma, direction = "long", idvar = "Subject",
        varying = matrix(colnames(plasma)[1:8], nr = 1),
        timevar = "time")
colnames(plasma)[4] <- "plasma"

toLatex(HSAURtable(plasma), pcol = 2,
    caption = "Plasma inorganic phosphate levels from 33 subjects, with data in ``long'' form.",
    label = "ch:LME:plasma:tab", rownames = FALSE)  


###################################################
### code chunk number 8: ch:LME:timber:lm
###################################################
summary(lm(loads ~ slippage, data = timber))


###################################################
### code chunk number 9: ch:LME:timber:plot (eval = FALSE)
###################################################
## xyplot(loads ~ slippage | specimen, data = timber,
##        layout = c(4, 2))


###################################################
### code chunk number 10: ch:LME:timber:plot
###################################################
plot(xyplot(loads ~ slippage | specimen, data = timber,
            layout = c(4, 2)))


###################################################
### code chunk number 11: ch:LME:timber:lme
###################################################
timber.lme <- lme(loads ~ slippage, 
                  random = ~1 | specimen,
                  data = timber, method = "ML")
timber.lme1 <- lme(loads ~ slippage, 
                   random = ~slippage | specimen,
                   data = timber, method = "ML")


###################################################
### code chunk number 12: ch:LME:timber:LRT
###################################################
library("RLRsim")
exactRLRT(timber.lme)


###################################################
### code chunk number 13: ch:LME:timber:anova
###################################################
anova(timber.lme, timber.lme1)


###################################################
### code chunk number 14: ch:LME:timber:lme:summary
###################################################
summary(timber.lme1)


###################################################
### code chunk number 15: ch:LME:timber:lme:predict
###################################################
timber$pred1 <- predict(timber.lme1)


###################################################
### code chunk number 16: ch:LME:timber:plotpredict (eval = FALSE)
###################################################
## pfun <- function(x, y) {
##     panel.xyplot(x, y[1:length(x)])
##     panel.lines(x, y[1:length(x) + length(x)], lty = 1)
## }
## plot(xyplot(cbind(loads, pred1) ~ slippage | specimen, 
##             data = timber, panel = pfun, layout = c(4, 2), 
##             ylab = "loads"))


###################################################
### code chunk number 17: ch:LME:timber:plotpredict
###################################################
pfun <- function(x, y) {
    panel.xyplot(x, y[1:length(x)])
    panel.lines(x, y[1:length(x) + length(x)], lty = 1)
}
plot(xyplot(cbind(loads, pred1) ~ slippage | specimen, data = timber,
       panel = pfun, layout = c(4, 2), ylab = "loads"))


###################################################
### code chunk number 18: ch:LME:timber:quad
###################################################
timber.lme2 <- lme(loads ~ slippage + I(slippage^2),
                   random = ~slippage | specimen,
                   data = timber, method = "ML")
anova(timber.lme1, timber.lme2)


###################################################
### code chunk number 19: ch:LME:timber:quad:plot
###################################################
timber$pred2 <- predict(timber.lme2)
plot(xyplot(cbind(loads, pred2) ~ slippage | specimen, 
            data = timber, panel = pfun, layout = c(4, 2), 
            ylab = "loads"))


###################################################
### code chunk number 20: ch:LME:plasma:plot (eval = FALSE)
###################################################
## x <- reshape(plasma, direction = "wide", timevar = "time", 
##              idvar = "Subject", v.names = "plasma")
## parallel(~ x[,-(1:2)] | group, data = x, horizontal = FALSE,
##          col = "black", scales = list(x = list(labels = 1:8)),
##          ylab = "Plasma inorganic phosphate",
##          xlab = "Time (hours after oral glucose challenge)")


###################################################
### code chunk number 21: ch:LME:plasma:plot
###################################################
x <- reshape(plasma, direction = "wide", timevar = "time", 
             idvar = "Subject", v.names = "plasma")
plot(parallel(~ x[,-(1:2)] | group, data = x, horizontal = FALSE,
         col = "black", scales = list(x = list(labels = 1:8)),
         ylab = "Plasma inorganic phosphate",
         xlab = "Time (hours after oral glucose challenge)"))


###################################################
### code chunk number 22: ch:LME:plasma:splom
###################################################
plot(splom(~ x[, grep("plasma", colnames(x))] | group, data = x, 
           cex = 1.5, pch = ".", pscales = NULL, varnames = 1:8))


###################################################
### code chunk number 23: ch:LME:plasma:lme
###################################################
plasma.lme1 <- lme(plasma ~ time + I(time^2) + group,
                   random = ~ time | Subject,   
                   data = plasma, method = "ML")
summary(plasma.lme1)


###################################################
### code chunk number 24: ch:LME:plasma:lm
###################################################
summary(lm(plasma ~ time + I(time^2) + group, data = plasma))


###################################################
### code chunk number 25: ch:LME:plasma:predictplot1
###################################################
pfun <- function(x, y, subscripts, groups) {  
    panel.xyplot(x, y[1:length(x)], pch = c(1:2)[groups[subscripts]])
    panel.lines(x, y[1:length(x) + length(x)], lty = 1)
}
plasma$pred1 <- predict(plasma.lme1)
plot(xyplot(cbind(plasma, pred1) ~ time | Subject, data = plasma, groups = group, 
       type = "b",
       pch = c(1, 2), ylab = "Plasma inorganic phosphate",
       xlab = "Time (hours after oral glucose challenge)", panel = pfun))


###################################################
### code chunk number 26: ch:LME:plasma:inter
###################################################
plasma.lme2 <- lme(plasma ~ time*group +I(time^2),
                   random = ~time | Subject, 
                   data = plasma, method = "ML")


###################################################
### code chunk number 27: ch:LME:plasma:interaov
###################################################
anova(plasma.lme1, plasma.lme2) 


###################################################
### code chunk number 28: ch:LME:plasma:summary
###################################################
summary(plasma.lme2)


###################################################
### code chunk number 29: ch:LME:plasma:predictplot2
###################################################
plasma$pred2 <- predict(plasma.lme2)
plot(xyplot(cbind(plasma, pred2) ~ time | Subject, data = plasma, groups = group,
       type = "b",
       pch = c(1, 2), ylab = "Plasma inorganic phosphate",
       xlab = "Time (hours after oral glucose challenge)", panel = pfun))


###################################################
### code chunk number 30: ch:LME:plasma:random
###################################################
res.int <- random.effects(plasma.lme2)[,1]
res.slope <- random.effects(plasma.lme2)[,2]


###################################################
### code chunk number 31: ch:LME:plasma:fitplot
###################################################
qqnorm(res.int,ylab="Estimated random intercepts",main="Random intercepts")
qqnorm(res.slope,ylab="Estimated random slopes",main="Random slopes")
resids<-resid(plasma.lme2) 
qqnorm(resids,ylab="Estimated residuals",main="Residuals")


###################################################
### code chunk number 32: ch:LME:BtheB:tab
###################################################
data("BtheB", package = "HSAUR2")
toLatex(HSAURtable(BtheB), pcol = 1,
    caption = "Data of a randomised trial evaluating the effects of Beat the Blues.",
    label = "ch:LME:BtheB:tab")


###################################################
### code chunk number 33: ch:LME:BtheB:plot
###################################################
ylim <- range(BtheB[,grep("bdi", names(BtheB))],
              na.rm = TRUE)
tau <- subset(BtheB, treatment == "TAU")[,
    grep("bdi", names(BtheB))]
boxplot(tau, main = "Treated as Usual", ylab = "BDI",
        xlab = "Time (in months)", names = c(0, 2, 3, 5, 8),
        ylim = ylim)
btheb <- subset(BtheB, treatment == "BtheB")[,
    grep("bdi", names(BtheB))]
boxplot(btheb, main = "Beat the Blues", ylab = "BDI",
        xlab = "Time (in months)", names = c(0, 2, 3, 5, 8),
        ylim = ylim)


###################################################
### code chunk number 34: ch:LME:BtheB:long
###################################################
BtheB$subject <- factor(rownames(BtheB))
nobs <- nrow(BtheB)
BtheB_long <- reshape(BtheB, idvar = "subject",
    varying = c("bdi.2m", "bdi.3m", "bdi.5m", "bdi.8m"),
    direction = "long")
BtheB_long$time <- rep(c(2, 3, 5, 8), rep(nobs, 4))


###################################################
### code chunk number 35: ch:LME:BtheB:lme
###################################################
BtheB_lme1 <- lme(bdi ~ bdi.pre + time + treatment + drug +
    length, random = ~ 1 | subject, data = BtheB_long, 
    na.action = na.omit)
BtheB_lme2 <- lme(bdi ~ bdi.pre + time + treatment + drug +
    length, random = ~ time | subject, data = BtheB_long,  
    na.action = na.omit)   


###################################################
### code chunk number 36: ch:LME:BtheB:aov
###################################################
anova(BtheB_lme1, BtheB_lme2)


###################################################
### code chunk number 37: ch:LME:BtheB:summary
###################################################
summary(BtheB_lme1)


###################################################
### code chunk number 38: ch:LME:BtheB-missing
###################################################
bdi <- BtheB[, grep("bdi", names(BtheB))]
plot(1:4, rep(-0.5, 4), type = "n", axes = FALSE, 
     ylim = c(0, 50), xlab = "Months", ylab = "BDI")
axis(1, at = 1:4, labels = c(0, 2, 3, 5))
axis(2)
for (i in 1:4) {
    dropout <- is.na(bdi[,i + 1])
    points(rep(i, nrow(bdi)) + ifelse(dropout, 0.05, -0.05), 
           jitter(bdi[,i]), pch = ifelse(dropout, 20, 1))
}


