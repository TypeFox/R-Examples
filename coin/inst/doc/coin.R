### R code from vignette source 'coin.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width = 60)
require("coin")
set.seed(290875)


###################################################
### code chunk number 2: YOY-kruskal
###################################################
library("coin")
YOY <- data.frame(
    length = c(46, 28, 46, 37, 32, 41, 42, 45, 38, 44,
               42, 60, 32, 42, 45, 58, 27, 51, 42, 52,
               38, 33, 26, 25, 28, 28, 26, 27, 27, 27,
               31, 30, 27, 29, 30, 25, 25, 24, 27, 30),
    site = gl(4, 10, labels = as.roman(1:4))
)

it <- independence_test(length ~ site, data = YOY,
          ytrafo = function(data)
              trafo(data, numeric_trafo = rank_trafo),
          teststat = "quadratic")
it


###################################################
### code chunk number 3: YOY-T
###################################################
statistic(it, type = "linear")


###################################################
### code chunk number 4: YOY-EV
###################################################
expectation(it)
covariance(it)


###################################################
### code chunk number 5: YOY-S
###################################################
statistic(it, type = "standardized")


###################################################
### code chunk number 6: YOY-c
###################################################
statistic(it)


###################################################
### code chunk number 7: YOY-p
###################################################
pvalue(it)


###################################################
### code chunk number 8: YOY-KW
###################################################
kt <- kruskal_test(length ~ site, data = YOY,
          distribution = approximate(B = 10000))
kt


###################################################
### code chunk number 9: YOY-KWp
###################################################
pvalue(kt)


###################################################
### code chunk number 10: jobsatisfaction-cmh
###################################################
data("jobsatisfaction", package = "coin")

ct <- cmh_test(jobsatisfaction)
ct


###################################################
### code chunk number 11: jobsatisfaction-s
###################################################
statistic(ct, type = "standardized")


###################################################
### code chunk number 12: jobsatisfaction-lbl
###################################################
lbl_test(jobsatisfaction)


###################################################
### code chunk number 13: jobsatisfaction-lbl-sc
###################################################
lbl_test(jobsatisfaction,
    scores = list(Job.Satisfaction = c(1, 3, 4, 5),
                  Income = c(3, 10, 20, 35)))


###################################################
### code chunk number 14: jobsatisfaction-lbl-sc-alt
###################################################
lbl_test(jobsatisfaction,
    ytrafo = function(data)
        trafo(data, ordered_trafo = function(y)
            of_trafo(y, scores = c(1, 3, 4, 5))),
    xtrafo = function(data)
        trafo(data, ordered_trafo = function(x)
            of_trafo(x, scores = c(3, 10, 20, 35))))


###################################################
### code chunk number 15: eggs-Durbin
###################################################
egg_data <- data.frame(
    scores = c(9.7, 8.7, 5.4, 5.0, 9.6, 8.8, 5.6,  3.6, 9.0,
               7.3, 3.8, 4.3, 9.3, 8.7, 6.8, 3.8, 10.0, 7.5,
               4.2, 2.8, 9.6, 5.1, 4.6, 3.6, 9.8,  7.4, 4.4,
               3.8, 9.4, 6.3, 5.1, 2.0, 9.4, 9.3,  8.2, 3.3,
               8.7, 9.0, 6.0, 3.3, 9.7, 6.7, 6.6,  2.8, 9.3,
               8.1, 3.7, 2.6, 9.8, 7.3, 5.4, 4.0,  9.0, 8.3,
               4.8, 3.8, 9.3, 8.3, 6.3, 3.8),
    sitting = factor(rep(c(1:15), rep(4, 15))),
    product = factor(c(1, 2, 4,  5, 2, 3, 6, 10, 2, 4, 6,  7,
                       1, 3, 5,  7, 1, 4, 8, 10, 2, 7, 8,  9,
                       2, 5, 8, 10, 5, 7, 9, 10, 1, 2, 3,  9,
                       4, 5, 6,  9, 1, 6, 7, 10, 3, 4, 9, 10,
                       1, 6, 8,  9, 3, 4, 7,  8, 3, 5, 6,  8))
)

independence_test(scores ~ product | sitting,
    data = egg_data, teststat = "quadratic",
    ytrafo = function(data)
        trafo(data, numeric_trafo = rank_trafo,
              block = egg_data$sitting))


###################################################
### code chunk number 16: eggs-Durbin-approx
###################################################
pvalue(independence_test(scores ~ product | sitting,
           data = egg_data, teststat = "quadratic",
           ytrafo = function(data)
               trafo(data, numeric_trafo = rank_trafo,
                     block = egg_data$sitting),
           distribution = approximate(B = 19999)))


###################################################
### code chunk number 17: eggs-Page
###################################################
independence_test(scores ~ product | sitting,
    data = egg_data,
    ytrafo = function(data)
        trafo(data, numeric_trafo = rank_trafo,
              block = egg_data$sitting),
    scores = list(product = 1:10))


###################################################
### code chunk number 18: warpbreaks-Tukey
###################################################
it <- independence_test(length ~ site, data = YOY,
          xtrafo = mcp_trafo(site = "Tukey"),
          teststat = "maximum",
          distribution = "approximate")
pvalue(it)
pvalue(it, method = "single-step")


