library(RTDE)

# ?EPD

#####
# (1) density function
x <- seq(0, 5, length=24)

cbind(x, dEPD(x, 1/2, 1/4, -1))

#####
# (2) distribution function

cbind(x, pEPD(x, 1/2, 1/4, -1, lower=FALSE))
cbind(x, pEPD(x, 1/2, 1/4, -1))

#####
# (3) quantile function

qEPD(1/2, 1/2, 1/4, -1)
qEPD(1:10/11, 1/2, 1/4, -1)
cbind(x, qEPD(pEPD(x, 1/2, 1/4, -1), 1/2, 1/4, -1)) #first five lines should be different


x <- seq(1, 10, length=21)
cbind(x, qEPD(pEPD(x, 1/2, 1/4, -1), 1/2, 1/4, -1))

system.time(x <- qEPD(1:10000/10001, 1/2, 1/4, -1))
head(x)
tail(x)

