require(lokern)
data(xSim)
(n <- length(xSim))
summary(xSim)
tt <- ((1:n) - 1/2)/n # equidistant x

str(lk <- lokerns(tt, xSim, trace=TRUE))
summary(lk$est)
summary((lk1 <- lokerns(tt,xSim, deriv = 1))$est)
summary((lk2 <- lokerns(tt,xSim, deriv = 2, trace=TRUE))$est)

summary(lk  $bandwidth)
summary(lk1 $bandwidth)
summary(lk2 $bandwidth)

str(lkH <- lokerns(tt, xSim, hetero=TRUE))
