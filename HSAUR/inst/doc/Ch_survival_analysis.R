### R code from vignette source 'Ch_survival_analysis.Rnw'
### Encoding: ASCII

###################################################
### code chunk number 1: setup
###################################################
rm(list = ls())
s <- search()[-1]
s <- s[-match(c("package:base", "package:stats", "package:graphics", "package:grDevices",
                "package:utils", "package:datasets", "package:methods", "Autoloads"), s)]
if (length(s) > 0) sapply(s, detach, character.only = TRUE)
if (!file.exists("tables")) dir.create("tables")
if (!file.exists("figures")) dir.create("figures")
set.seed(290875)
options(prompt = "R> ", continue = "+  ",
    width = 63, # digits = 4, 
    SweaveHooks = list(leftpar = function() 
        par(mai = par("mai") * c(1, 1.05, 1, 1))))
HSAURpkg <- require("HSAUR")
if (!HSAURpkg) stop("cannot load package ", sQuote("HSAUR"))
rm(HSAURpkg)
a <- Sys.setlocale("LC_ALL", "C")
book <- TRUE
refs <- cbind(c("AItR", "SI", "CI", "ANOVA", "MLR", "GLM", 
                "DE", "RP", "SA", "ALDI", "ALDII", "MA", "PCA", 
                "MDS", "CA"), 1:15)
ch <- function(x, book = TRUE) {
    ch <- refs[which(refs[,1] == x),]
    if (book) {
        return(paste("Chapter~\\\\ref{", ch[1], "}", sep = ""))
    } else {
        return(paste("Chapter~\\\\ref{", ch[2], "}", sep = ""))
    }
}


###################################################
### code chunk number 2: SA-setup
###################################################
x <- library("survival")
x <- library("coin")
x <- library("party")


###################################################
### code chunk number 3: SA-glioma-KM
###################################################
data("glioma", package = "coin")
library("survival")
layout(matrix(1:2, ncol = 2))
g3 <- subset(glioma, histology == "Grade3")
plot(survfit(Surv(time, event) ~ group, data = g3), 
     main = "Grade III Glioma", lty = c(2, 1), 
     ylab = "Probability", xlab = "Survival Time in Month",
     legend.bty = "n", legend.text = c("Control", "Treated")
)
g4 <- subset(glioma, histology == "GBM")
plot(survfit(Surv(time, event) ~ group, data = g4), 
     main = "Grade IV Glioma", ylab = "Probability", 
     lty = c(2, 1), xlab = "Survival Time in Month", 
     xlim = c(0, max(glioma$time) * 1.05))


###################################################
### code chunk number 4: SA-glioma-logrank
###################################################
survdiff(Surv(time, event) ~ group, data = g3)


###################################################
### code chunk number 5: SA-glioma-exact
###################################################
library("coin")
surv_test(Surv(time, event) ~ group, data = g3, 
          distribution = "exact")


###################################################
### code chunk number 6: SA-glioma-g4
###################################################
surv_test(Surv(time, event) ~ group, data = g4, 
                  distribution = "exact")


###################################################
### code chunk number 7: SA-glioma-hist
###################################################
surv_test(Surv(time, event) ~ group | histology, data = glioma,
          distribution = approximate(B = 10000))


###################################################
### code chunk number 8: SA-GBSG2-plot
###################################################
data("GBSG2", package = "TH.data")
plot(survfit(Surv(time, cens) ~ horTh, data = GBSG2), 
     lty = 1:2, mark.time = FALSE,  ylab = "Probability", 
     xlab = "Survival Time in Days")
legend(250, 0.2, legend = c("yes", "no"), lty = c(2, 1), 
       title = "Hormonal Therapy", bty = "n")


###################################################
### code chunk number 9: SA-GBSG2-coxph
###################################################
GBSG2_coxph <- coxph(Surv(time, cens) ~ ., data = GBSG2)


###################################################
### code chunk number 10: SA-GBSG2-coxph-ci
###################################################
ci <- confint(GBSG2_coxph)
exp(cbind(coef(GBSG2_coxph), ci))["horThyes",]


###################################################
### code chunk number 11: GBSG2-coxph-summary
###################################################
summary(GBSG2_coxph)


###################################################
### code chunk number 12: SA-GBSG2-zph
###################################################
GBSG2_zph <- cox.zph(GBSG2_coxph)
GBSG2_zph


###################################################
### code chunk number 13: SA-GBSG2-zph-plot
###################################################
plot(GBSG2_zph, var = "age")


###################################################
### code chunk number 14: SA-GBSG2-Martingal
###################################################
layout(matrix(1:3, ncol = 3))
res <- residuals(GBSG2_coxph)
plot(res ~ age, data = GBSG2, ylim = c(-2.5, 1.5), 
     pch = ".", ylab = "Martingale Residuals")
abline(h = 0, lty = 3)
plot(res ~ pnodes, data = GBSG2, ylim = c(-2.5, 1.5), 
     pch = ".", ylab = "")
abline(h = 0, lty = 3)
plot(res ~ log(progrec), data = GBSG2, ylim = c(-2.5, 1.5), 
     pch = ".", ylab = "")
abline(h = 0, lty = 3)


###################################################
### code chunk number 15: SA-GBSG2-ctree
###################################################
GBSG2_ctree <- ctree(Surv(time, cens) ~ ., data = GBSG2)


###################################################
### code chunk number 16: SA-GBSG2-ctree-plot
###################################################
plot(GBSG2_ctree)


