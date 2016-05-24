### R code from vignette source 'BGPhazard.Rnw'

###################################################
### code chunk number 1: BGPhazard.Rnw:235-238 (eval = FALSE)
###################################################
## data(gehan)
## times <- gehan$time[gehan$treat == "6-MP"]
## delta <- gehan$cens[gehan$treat == "6-MP"]


###################################################
### code chunk number 2: BGPhazard.Rnw:247-250 (eval = FALSE)
###################################################
## ExG1 <- GaMRes(times, delta, type.t = 2, K = 35, type.c = 1, 
##                iterations = 3000)
## GaPloth(ExG1, confint = FALSE)


###################################################
### code chunk number 3: BGPhazard.Rnw:270-273 (eval = FALSE)
###################################################
## ExG2 <- GaMRes(times, delta, type.t = 2, K = 35, type.c = 2, 
##                c.r = rep(50, 34), iterations = 3000)
## GaPloth(ExG2)


###################################################
### code chunk number 4: BGPhazard.Rnw:291-292 (eval = FALSE)
###################################################
## GaPlotDiag(ExG2, variable = "lambda", pos = 6)


###################################################
### code chunk number 5: BGPhazard.Rnw:308-311 (eval = FALSE)
###################################################
## ExG3 <- GaMRes(times, delta, type.t = 2, K = 35, type.c = 3, epsilon = 1,
##                iterations = 3000)
## GaPloth(ExG3)


###################################################
### code chunk number 6: BGPhazard.Rnw:332-335 (eval = FALSE)
###################################################
## ExG4 <- GaMRes(times, delta, type.t = 2, K = 35, type.c = 4, 
##                iterations = 3000)
## GaPloth(ExG4)


###################################################
### code chunk number 7: BGPhazard.Rnw:356-359 (eval = FALSE)
###################################################
## ExG5 <- GaMRes(times, delta, type.t = 1, K = 8, type.c = 2, 
##                c.r = rep(50, 7), iterations = 3000)
## GaPloth(ExG5)


###################################################
### code chunk number 8: BGPhazard.Rnw:380-383 (eval = FALSE)
###################################################
## ExG6 <- GaMRes(times, delta, type.t = 1, K = 8, type.c = 3,
##                iterations=3000)
## GaPloth(ExG6)


###################################################
### code chunk number 9: BGPhazard.Rnw:404-407 (eval = FALSE)
###################################################
## ExG7 <- GaMRes(times, delta, type.t = 1, K = 8, type.c = 4,
##                iterations = 3000)
## GaPloth(ExG7)


###################################################
### code chunk number 10: BGPhazard.Rnw:429-432 (eval = FALSE)
###################################################
## data(psych)
## times <- psych$time
## delta <- psych$death


###################################################
### code chunk number 11: BGPhazard.Rnw:439-441 (eval = FALSE)
###################################################
## ExB1 <- BeMRes(times, delta, type.c = 1, iterations = 3000)
## BePloth(ExB1, confint = FALSE)


###################################################
### code chunk number 12: BGPhazard.Rnw:462-465 (eval = FALSE)
###################################################
## ExB2 <- BeMRes(times, delta, type.c = 2, c.r = rep(100, 39), 
##                iterations = 3000)
## BePloth(ExB2)


###################################################
### code chunk number 13: BGPhazard.Rnw:484-485 (eval = FALSE)
###################################################
## BePlotDiag(ExB2, variable = "Pi", pos = 6)


###################################################
### code chunk number 14: BGPhazard.Rnw:499-501 (eval = FALSE)
###################################################
## ExB3 <- BeMRes(times, delta, type.c = 3, epsilon = 1, iterations = 3000)
## BePloth(ExB3)


###################################################
### code chunk number 15: BGPhazard.Rnw:522-524 (eval = FALSE)
###################################################
## ExB4 <- BeMRes(times, delta, type.c = 4, iterations = 3000)
## BePloth(ExB4)


###################################################
### code chunk number 16: BGPhazard.Rnw:579-601 (eval = FALSE)
###################################################
## SampWeibull <- function(n, a = 10, b = 1, beta = c(1, 1)) {
##   M <- matrix(0, ncol = 7, nrow = n)
##   for(i in 1:n){
##     M[i, 1] <- i
##     M[i, 2] <- x1 <- runif(1)
##     M[i, 3] <- x2 <- runif(1)
##     M[i, 4] <- rweibull(1, shape = b, 
##                         scale = 1 / (a * exp(cbind(x1, x2) %*% beta)))
##     M[i, 5] <- rexp(1)
##     M[i, 6] <- M[i, 4] > M[i, 5]
##     M[i, 7] <- min(M[i, 4], M[i, 5])
##   }
##   colnames(M) <- c("i", "x_i1", "x_i2", "t_i", "c_i", "delta",
##                    "min{c_i, d_i}")
##   return(M)
## }
## dat <- SampWeibull(100, 0.1, 1, c(1, 1))
## dat <- cbind(dat[, c(4, 6)], dat[, c(2, 3)])
## CG <- CGaMRes(dat, K = 10, iterations = 3000, thpar = 10)
## CGaPloth(CG)
## PlotTheta(CG)
## CGaPred(CG)


