library(RTDE)

# ?EPD

#####
# (1) density function
x <- seq(0, 5, length=24)

cbind(x, dfrechet(x, 1/2, 1/4))

cbind(x, dfrechet(x, 1, 0))

#####
# (2) distribution function

cbind(x, pfrechet(x, 1/2, 1/4))
cbind(x, pfrechet(x, 1, 0))

#####
# (3) quantile function

p <- 1:10/11
cbind(p, qfrechet(p, 1, 0))
