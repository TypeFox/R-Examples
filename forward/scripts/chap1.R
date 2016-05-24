#-----------------------------------------------------------------------------#
#   Chapter 1 
#-----------------------------------------------------------------------------#

data(forbes)
plot(forbes, xlab="Boiling point", ylab="100 × log(pressure)")
mod <- fwdlm(y ~ x, data=forbes)
summary(mod)
plot(mod)

#-----------------------------------------------------------------------------#

data(ar)
pairs(ar)
mod <- fwdlm(y ~ x1 + x2 + x3, data=ar)
summary(mod)
plot(mod)
plot(mod, squared=T, 1) 

#-----------------------------------------------------------------------------#

data(wool)

mod.ols <- lm(y ~ x1 + x2 + x3, data=wool)
plot(mod.ols)
qqnorm(mod.ols$residuals)

library(MASS)
boxcox(mod.ols, plotit=T)

mod <- fwdsco(y ~ x1 + x2 + x3, lambda=0, data=wool)
summary(mod)
plot(mod, ylim=c(-10,10), plot.mle=FALSE)

mod <- fwdlm(y ~ x1 + x2 + x3, data=wool)
summary(mod)
plot(mod) 
plot(mod, 1)
plot(mod, 10, ylim=c(0.7,1))
plot(mod, 5)

mod <- fwdlm(log(y) ~ x1 + x2 + x3, data=wool)
summary(mod)
plot(mod) 
plot(mod, 1, ylim=c(-4.5, 2))
plot(mod, 5, scaled=FALSE, ylim=c(-2,6.5))
