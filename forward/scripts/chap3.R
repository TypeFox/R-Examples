#-----------------------------------------------------------------------------#
#   Chapter 3 
#-----------------------------------------------------------------------------#

data("hawkins")
pairs(hawkins)
mod.ols <- lm(y ~ ., data=hawkins)
qqnorm(mod.ols$residuals, ylab="LS residuals")

mod.lms <- lmsreg(y ~ ., data=hawkins)
qqnorm(mod.lms$residuals, ylab="LMS residuals")

mod <- fwdlm(y ~ ., data=hawkins)
summary(mod)
plot(mod) 
plot(mod, 1, squared=TRUE, xlim=c(0,132))

par(mfrow=c(1,2))
plot(mod, 5, scaled=FALSE)
plot(mod, 6, ylim=c(-40,40))

par(mfrow=c(1,1))
plot(mod, 8, ylim=c(0,20))

par(mfrow=c(2,1))
plot(mod, 3, ylim=c(0,10))
plot(mod, 4, ylim=c(0,20))

par(mfrow=c(1,1))
plot(mod, 9)

n <- nrow(hawkins)
inc <- mod$included
pch <- rep(1,n)
pch[inc$"m=86"] <- 3                         # first 86 obs to enter
pch[setdiff(inc$"m=128", inc$"m=122")] <- 15 # last 6 obs to enter
pch[setdiff(inc$"m=122", inc$"m=110")] <- 17 # obs which enter between step 110 and 122
par(mfrow=c(1,1))
plot(y ~ x8, data=hawkins, pch=pch)

pairs(hawkins, col=as.numeric(factor(pch)), pch=pch)

#-----------------------------------------------------------------------------#

data(stackloss)
pairs(stackloss)
mod1.ols <- lm(Loss ~ . ,data=stackloss)
summary(mod1.ols)
plot(mod1.ols) 

mod1 <- fwdlm(Loss ~ ., data=stackloss, nsamp="exact")
summary(mod1)
plot(mod1) 
plot(mod1, 1, squared=T)
plot(mod1, 2)
plot(mod1, 6, ylim=c(-40,40))

mod2 <- fwdlm(Loss ~ Air + Temp, data=stackloss)
summary(mod2)
plot(mod2, 2)
plot(mod2, 10, ylim=c(0.9,1))

mod3.ols <- lm(Loss ~ Air*Temp + I(Air^2), data=stackloss)
summary(mod3.ols)

mod4 <- fwdlm(Loss ~ Air*Temp + I(Air^2), data=stackloss, nsamp="exact")
plot(mod4, 1) 
plot(mod4, 2)

mod5 <- fwdsco(Loss ~ Air*Temp + I(Air^2), data=stackloss, lambda = 1, nsamp="exact")
summary(mod5)
plot(mod5)

mod6 <- fwdlm(log(Loss) ~ Air*Temp + I(Air^2), data=stackloss, nsamp="exact")
summary(mod6)
plot(mod6, 1)
plot(mod6, 6)
plot(mod6, 8)

mod7 <- fwdlm(log(Loss) ~ Air + Temp, data=stackloss, nsamp="exact")
plot(mod7, 10, ylim=c(0.9,1))
plot(mod7, 1, ylim=c(-4,4))

mod8 <- fwdsco(Loss ~ Air + Temp, data=stackloss, lambda = 0, nsamp="exact")
summary(mod8)
plot(mod8, ylim=c(-4,4))

mod9 <- fwdsco(Loss ~ Air + Temp, data=stackloss, lambda = 0.5, nsamp="exact")
plot(mod9, ylim=c(-4,4))

mod10 <- fwdlm(sqrt(Loss) ~ Air + Temp, data=stackloss, nsamp="exact")
summary(mod10)
plot(mod10, 6, ylim=c(-40,40))
plot(mod10, 10, ylim=c(0.9,1))
plot(mod10, 1, ylim=c(-4,4))


#-----------------------------------------------------------------------------#

data(salinity)

mod1 <- fwdlm(y ~ x1 + x2 + x3, data=salinity, nsamp="exact")
plot(mod1)
plot(mod1, 1)
plot(mod1, 2)

mark <- rep(16, nrow(salinity)); mark[16] <- 1
pairs(salinity, pch=mark)

salinity$x3[16] <- 23.443  # correction asssuming misprint
mod2 <- fwdlm(y ~ x1 + x2 + x3, data=salinity, nsamp="exact")
plot(mod2, 1, th.Res=1.5)

par(mfrow=c(1,2))
plot(mod2, 6, ylim=c(-40,40))
mark <- rep(1, nrow(salinity)); mark[c(9, 15, 17)] <- 16
plot(y ~ x1, data=salinity, pch=mark, ylab="Salinity", xlab="Lagged salinity")
text(salinity$x1[mark==16], salinity$y[mark==16], which(mark==16), pos=4)
par(mfrow=c(1,1))

mod3 <- fwdlm(y ~ x1 + x3, data=salinity, nsamp="exact")
plot(mod3, 1)
plot(mod3, 6, ylim=c(-40,40))

col <- pch <- rep(1,nrow(salinity))
col[c(6,14,15)]  <- pch[c(6,14,15)]  <- 3  # my best start subset
col[c(9,15,17)] <- pch[c(9,15,17)] <- 4  # AR outliers
pairs(salinity[,c(4,1,3)], col=col, pch=pch, 
      upper.panel = function(x,y,...) text(x,y,1:length(x), ...) )

#-----------------------------------------------------------------------------#

data(ozone)

mod1 <- fwdlm(y ~ ., data=ozone)
plot(mod1, 1)
plot(mod1, 6,  ylim=c(-10,10))
plot(mod1, 8, ylim=c(0,6))

mod1.sco <- fwdsco(y ~ ., data=ozone, lambda=1)
plot(mod1.sco)

mod2 <- lm(log(y) ~ ., data=ozone)
Time <- 1:nrow(ozone)
library(MASS)
plot(Time, studres(mod2), type="l")
abline(0, 0, lty=2)

plot(fwdsco(y ~ ., data=ozone, lambda=0), ylim=c(-6,6))

ozone$Time <- Time
mod3 <- lm(log(y) ~ ., data=ozone)
summary(mod3)

mod4 <- lm(log(y) ~ Time + x2 + x5 + x6 + x8, data=ozone)
summary(mod4)

mod5 <- fwdlm(log(y) ~ Time + x2 + x5 + x6 + x8, data=ozone)
plot(mod5, 1, ylim=c(-4.5, 2.5))
plot(mod5, 2)

mod5.sco <- fwdsco(y ~ Time + x2 + x5 + x6 + x8, data=ozone, lambda=0) 
par(mfrow=c(1,2))
plot(mod5.sco, ylim=c(-6,6), plot.mle=F)
plot(mod5, 6, ylim=c(-40,40))

