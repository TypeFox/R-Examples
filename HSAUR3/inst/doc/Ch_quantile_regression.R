### R code from vignette source 'Ch_quantile_regression.Rnw'
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
    show.signif.stars = FALSE,
    SweaveHooks = list(leftpar = function() 
        par(mai = par("mai") * c(1, 1.05, 1, 1)),
        bigleftpar = function()
        par(mai = par("mai") * c(1, 1.7, 1, 1))))
HSAURpkg <- require("HSAUR3")
if (!HSAURpkg) stop("cannot load package ", sQuote("HSAUR2"))
rm(HSAURpkg)
 ### </FIXME> hm, R-2.4.0 --vanilla seems to need this
a <- Sys.setlocale("LC_ALL", "C")
 ### </FIXME>
book <- TRUE
refs <- cbind(c("AItR", "DAGD", "SI", "CI", "ANOVA", "MLR", "GLM", 
                "DE", "RP", "GAM", "SA", "ALDI", "ALDII", "SIMC", "MA", "PCA", 
                "MDS", "CA"), 1:18)
ch <- function(x) {
    ch <- refs[which(refs[,1] == x),]
    if (book) {
        return(paste("Chapter~\\\\ref{", ch[1], "}", sep = ""))
    } else {
        return(paste("Chapter~", ch[2], sep = ""))
    }
}
if (file.exists("deparse.R"))
    source("deparse.R")

setHook(packageEvent("lattice", "attach"), function(...) {
    lattice.options(default.theme = 
        function()
            standard.theme("pdf", color = FALSE))
    })


###################################################
### code chunk number 2: singlebook
###################################################
book <- FALSE


###################################################
### code chunk number 3: QR-setup
###################################################
library("lattice")
trellis.par.set(list(plot.symbol = list(col=1,pch=20, cex=0.7),
                     box.rectangle = list(col=1),
                     plot.line = list(col = 1, lwd = 1),
                     box.umbrella = list(lty=1, col=1),
                     strip.background = list(col = "white")))
ltheme <- canonical.theme(color = FALSE)     ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme)
data("db", package = "gamlss.data")
nboys <- with(db, sum(age > 2))


###################################################
### code chunk number 4: QR-db
###################################################
summary(db)
db$cut <- cut(db$age, breaks = c(2, 9, 23), 
              labels = c("2-9 yrs", "9-23 yrs"))


###################################################
### code chunk number 5: QR-db-plot
###################################################
db$cut <- cut(db$age, breaks = c(2, 9, 23), 
              labels = c("2-9 yrs", "9-23 yrs"))
xyplot(head ~ age | cut, data = db, xlab = "Age (years)", 
       ylab = "Head circumference (cm)",
       scales = list(x = list(relation = "free")),
       layout = c(2, 1), pch = 19, 
       col = rgb(.1, .1, .1, .1))


###################################################
### code chunk number 6: QR-db-lm2.9.23
###################################################
(lm2.9 <- lm(head ~ age, data = db, subset = age < 9))
(lm9.23 <- lm(head ~ age, data = db, subset = age > 9))


###################################################
### code chunk number 7: QR-db-lm
###################################################
(lm_mod <- lm(head ~ age:I(age < 9) + I(age < 9) - 1, 
              data = db))


###################################################
### code chunk number 8: QR-db-median
###################################################
library("quantreg")
(rq_med2.9 <- rq(head ~ age, data = db, tau = 0.5, 
                 subset = age < 9))
(rq_med9.23 <- rq(head ~ age, data = db, tau = 0.5, 
                  subset = age > 9))


###################################################
### code chunk number 9: QR-db-lmrq2.9
###################################################
cbind(coef(lm2.9)[1], confint(lm2.9, parm = "(Intercept)"))
cbind(coef(lm2.9)[2], confint(lm2.9, parm = "age"))
summary(rq_med2.9, se = "rank")


###################################################
### code chunk number 10: QR-db-lmrq9.23
###################################################
cbind(coef(lm9.23)[1], confint(lm9.23, parm = "(Intercept)"))
cbind(coef(lm9.23)[2], confint(lm9.23, parm = "age"))
summary(rq_med9.23, se = "rank")


###################################################
### code chunk number 11: QR-db-tau
###################################################
tau <- c(.01, .1, .25, .5, .75, .9, .99)


###################################################
### code chunk number 12: QR-db-age
###################################################
gage <- c(2:9, 9:23)
i <- 1:8


###################################################
### code chunk number 13: QR-db-lm-fit_05
###################################################
idf <- data.frame(age = gage[i])
p <- predict(lm2.9, newdata = idf, level = 0.5,
             interval = "prediction")
colnames(p) <- c("0.5", "0.25", "0.75")
p


###################################################
### code chunk number 14: QR-db-lm-fit
###################################################
p <- cbind(p, predict(lm2.9, newdata = idf, level = 0.8,
                      interval = "prediction")[,-1])
colnames(p)[4:5] <- c("0.1", "0.9")
p <- cbind(p, predict(lm2.9, newdata = idf, level = 0.98, 
                      interval = "prediction")[,-1])
colnames(p)[6:7] <- c("0.01", "0.99")
p2.9 <- p[, c("0.01", "0.1", "0.25", "0.5", 
              "0.75", "0.9", "0.99")]
idf <- data.frame(age = gage[-i])
p <- predict(lm9.23, newdata = idf, level = 0.5, 
             interval = "prediction")
colnames(p) <- c("0.5", "0.25", "0.75")
p <- cbind(p, predict(lm9.23, newdata = idf, level = 0.8, 
                      interval = "prediction")[,-1])
colnames(p)[4:5] <- c("0.1", "0.9")
p <- cbind(p, predict(lm9.23, newdata = idf, level = 0.98, 
                      interval = "prediction")[,-1])
colnames(p)[6:7] <- c("0.01", "0.99")


###################################################
### code chunk number 15: QR-db-lm-fit2
###################################################
p9.23 <- p[, c("0.01", "0.1", "0.25", "0.5", 
               "0.75", "0.9", "0.99")]
round((q2.23 <- rbind(p2.9, p9.23)), 3)


###################################################
### code chunk number 16: QR-db-lm-plot
###################################################
pfun <- function(x, y, ...) {
    panel.xyplot(x = x, y = y, ...)
    if (max(x) <= 9) {
        apply(q2.23, 2, function(x) 
              panel.lines(gage[i], x[i]))
    } else {
        apply(q2.23, 2, function(x) 
              panel.lines(gage[-i], x[-i]))
    }
    panel.text(rep(max(db$age), length(tau)), 
               q2.23[nrow(q2.23),], label = tau, cex = 0.9) 
    panel.text(rep(min(db$age), length(tau)), 
               q2.23[1,], label = tau, cex = 0.9)
}
xyplot(head ~ age | cut, data = db, xlab = "Age (years)", 
       ylab = "Head circumference (cm)", pch = 19,
       scales = list(x = list(relation = "free")),
       layout = c(2, 1), col = rgb(.1, .1, .1, .1), 
       panel = pfun)


###################################################
### code chunk number 17: QR-db-rq2.9
###################################################
(rq2.9 <- rq(head ~ age, data = db, tau = tau, 
             subset = age < 9))


###################################################
### code chunk number 18: QR-db-rq9.23
###################################################
(rq9.23 <- rq(head ~ age, data = db, tau = tau, 
              subset = age > 9))


###################################################
### code chunk number 19: QR-db-rq-fit
###################################################
p2.23 <- rbind(predict(rq2.9, 
                   newdata = data.frame(age = gage[i])),
               predict(rq9.23, 
                   newdata = data.frame(age = gage[-i])))


###################################################
### code chunk number 20: QR-db-rq-plot
###################################################
pfun <- function(x, y, ...) {
    panel.xyplot(x = x, y = y, ...)
    if (max(x) <= 9) {
        apply(q2.23, 2, function(x) 
              panel.lines(gage[i], x[i], lty = 2))
        apply(p2.23, 2, function(x) 
              panel.lines(gage[i], x[i]))
    } else {
        apply(q2.23, 2, function(x) 
              panel.lines(gage[-i], x[-i], lty = 2))
        apply(p2.23, 2, function(x) 
              panel.lines(gage[-i], x[-i]))
    }
    panel.text(rep(max(db$age), length(tau)), 
               p2.23[nrow(p2.23),], label = tau, cex = 0.9)
    panel.text(rep(min(db$age), length(tau)), 
               p2.23[1,], label = tau, cex = 0.9)
}
xyplot(head ~ age | cut, data = db, xlab = "Age (years)",
       ylab = "Head circumference (cm)", pch = 19,
       scales = list(x = list(relation = "free")),
       layout = c(2, 1), col = rgb(.1, .1, .1, .1),
       panel = pfun)


###################################################
### code chunk number 21: QR-db-rqss-fit
###################################################
rqssmod <- vector(mode = "list", length = length(tau))
db$lage <- with(db, age^(1/3))
for (i in 1:length(tau))
    rqssmod[[i]] <- rqss(head ~ qss(lage, lambda = 1),
                         data = db, tau = tau[i])


###################################################
### code chunk number 22: QR-db-rqss-pred
###################################################
gage <- seq(from = min(db$age), to = max(db$age), 
            length = 50)
p <- sapply(1:length(tau), function(i) {
    predict(rqssmod[[i]], 
        newdata = data.frame(lage = gage^(1/3)))
})


###################################################
### code chunk number 23: QR-db-rqss-plot
###################################################
pfun <- function(x, y, ...) {
    panel.xyplot(x = x, y = y, ...)
    apply(p, 2, function(x) panel.lines(gage, x))
    panel.text(rep(max(db$age), length(tau)), 
               p[nrow(p),], label = tau, cex = 0.9)
    panel.text(rep(min(db$age), length(tau)), 
               p[1,], label = tau, cex = 0.9)
}
xyplot(head ~ age | cut, data = db, xlab = "Age (years)",
       ylab = "Head circumference (cm)", pch = 19,
       scales = list(x = list(relation = "free")),
       layout = c(2, 1), col = rgb(.1, .1, .1, .1),
       panel = pfun)


