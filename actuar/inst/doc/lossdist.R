### R code from vignette source 'lossdist.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: lossdist.Rnw:53-55
###################################################
library(actuar)
options(width = 62, digits = 4)


###################################################
### code chunk number 2: lossdist.Rnw:265-268
###################################################
x <- grouped.data(Group = c(0, 25, 50, 100, 150, 250, 500),
                  Line.1 = c(30, 31, 57, 42, 65, 84),
                  Line.2 = c(26, 33, 31, 19, 16, 11))


###################################################
### code chunk number 3: lossdist.Rnw:272-273
###################################################
class(x)


###################################################
### code chunk number 4: lossdist.Rnw:277-278
###################################################
x


###################################################
### code chunk number 5: lossdist.Rnw:287-288
###################################################
x[, 1]                                  # group boundaries


###################################################
### code chunk number 6: lossdist.Rnw:292-293
###################################################
x[, -1]                                 # group frequencies


###################################################
### code chunk number 7: lossdist.Rnw:296-297
###################################################
x[1:3,]                                 # first 3 groups


###################################################
### code chunk number 8: lossdist.Rnw:306-308
###################################################
x[1, 2] <- 22; x                        # frequency replacement
x[1, c(2, 3)] <- c(22, 19); x           # frequency replacement


###################################################
### code chunk number 9: lossdist.Rnw:311-313
###################################################
x[1, 1] <- c(0, 20); x                  # boundary replacement
x[c(3, 4), 1] <- c(55, 110, 160); x


###################################################
### code chunk number 10: lossdist.Rnw:325-326
###################################################
mean(x)


###################################################
### code chunk number 11: lossdist.Rnw:336-337
###################################################
hist(x[, -3])


###################################################
### code chunk number 12: lossdist.Rnw:341-342
###################################################
hist(x[, -3])


###################################################
### code chunk number 13: lossdist.Rnw:377-381
###################################################
Fnt <- ogive(x)
knots(Fnt)                              # group boundaries
Fnt(knots(Fnt))                         # ogive at group boundaries
plot(Fnt)                               # plot of the ogive


###################################################
### code chunk number 14: lossdist.Rnw:385-386
###################################################
plot(Fnt)


###################################################
### code chunk number 15: lossdist.Rnw:395-397
###################################################
quantile(x)
Fnt(quantile(x))


###################################################
### code chunk number 16: lossdist.Rnw:408-410
###################################################
data("dental"); dental
data("gdental"); gdental


###################################################
### code chunk number 17: lossdist.Rnw:419-421
###################################################
emm(dental, order = 1:3)                # first three moments
emm(gdental, order = 1:3)               # idem


###################################################
### code chunk number 18: lossdist.Rnw:429-436
###################################################
lev <- elev(dental)
lev(knots(lev))                         # ELEV at data points
plot(lev, type = "o", pch = 19)         # plot of the ELEV function

lev <- elev(gdental)
lev(knots(lev))                         # ELEV at data points
plot(lev, type = "o", pch = 19)         # plot of the ELEV function


###################################################
### code chunk number 19: lossdist.Rnw:440-443
###################################################
par(mfrow = c(1, 2))
plot(elev(dental), type = "o", pch = 19)
plot(elev(gdental), type = "o", pch = 19)


###################################################
### code chunk number 20: lossdist.Rnw:510-511
###################################################
op <- options(warn = -1)                # hide warnings from mde()


###################################################
### code chunk number 21: lossdist.Rnw:513-516
###################################################
mde(gdental, pexp, start = list(rate = 1/200), measure = "CvM")
mde(gdental, pexp, start = list(rate = 1/200), measure = "chi-square")
mde(gdental, levexp, start = list(rate = 1/200), measure = "LAS")


###################################################
### code chunk number 22: lossdist.Rnw:518-519
###################################################
options(op)                             # restore warnings


###################################################
### code chunk number 23: lossdist.Rnw:525-527 (eval = FALSE)
###################################################
## mde(gdental, ppareto, start = list(shape = 3, scale = 600),
##         measure = "CvM") # no convergence


###################################################
### code chunk number 24: lossdist.Rnw:529-532
###################################################
out <- try(mde(gdental, ppareto, start = list(shape = 3, scale = 600),
        measure = "CvM"), silent = TRUE)
cat(sub(", measure", ",\n             measure", out))


###################################################
### code chunk number 25: lossdist.Rnw:538-542
###################################################
pparetolog <- function(x, logshape, logscale)
    ppareto(x, exp(logshape), exp(logscale))
( p <- mde(gdental, pparetolog, start = list(logshape = log(3),
                                logscale = log(600)), measure = "CvM") )


###################################################
### code chunk number 26: lossdist.Rnw:545-546
###################################################
exp(p$estimate)


###################################################
### code chunk number 27: lossdist.Rnw:642-648
###################################################
f <- coverage(pdf = dgamma, cdf = pgamma, deductible = 1, limit = 10)
f
f(0, shape = 5, rate = 1)
f(5, shape = 5, rate = 1)
f(9, shape = 5, rate = 1)
f(12, shape = 5, rate = 1)


###################################################
### code chunk number 28: lossdist.Rnw:659-662
###################################################
x <- rgamma(100, 2, 0.5)
y <- pmin(x[x > 1], 9)
op <- options(warn = -1)                # hide warnings from fitdistr()


###################################################
### code chunk number 29: lossdist.Rnw:664-666
###################################################
library(MASS)
fitdistr(y, f, start = list(shape = 2, rate = 0.5))


###################################################
### code chunk number 30: lossdist.Rnw:668-669
###################################################
options(op)                             # restore warnings


