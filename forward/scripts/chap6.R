#-----------------------------------------------------------------------------#
#   Chapter 6
#-----------------------------------------------------------------------------#

data(carinsuk)
carinsuk <- na.omit(carinsuk)

mod.glm <- glm(AvCost ~ OwnerAge + Model + CarAge, 
               weights = NClaims, data = carinsuk, 
               family=Gamma(inverse))
summary(mod.glm)
plot(mod.glm)

mod1 <- fwdglm(AvCost ~ OwnerAge + Model + CarAge, 
               weights = NClaims, data = carinsuk,
               family=Gamma(inverse), nsamp=10000)
# it takes a long time... and give as best subset
inibsb <- c(12, 21, 31, 49, 53, 55, 65, 77, 86, 89, 99, 101, 104, 109)
# ... so you may want to use
mod1 <- fwdglm(AvCost ~ OwnerAge + Model + CarAge, 
               weights=NClaims, data = carinsuk, 
               family=Gamma(inverse), bsb=inibsb)

summary(mod1, 5) # last five steps of forward-search

plot(mod1, 8)
plot(mod1, 1)
plot(mod1, 2)
plot(mod1, 7)
par(mfrow=c(1,2))
plot(mod1, 5)
plot(mod1, 6, ylim=c(-45,45))
par(mfrow=c(1,2))
plot(mod1, 10)
plot(mod1, 3)
par(mfrow=c(1,1))


#-----------------------------------------------------------------------------#

data(dialectric)

library(lattice)
trellis.device(bg="white")

xyplot(y ~ time | as.factor(temp), 
       data = dialectric, col=1,
       panel = function(x,y,...) 
               { panel.xyplot(x, y, type=c("p"), ...)
                 panel.loess(x, y, col=4)
                 panel.grid(lty=3, col="lightgrey")  })
xyplot(y ~ log(time) | as.factor(temp), 
       data = dialectric, col=1,
       panel = function(x,y,...) 
               { panel.xyplot(x, y, type=c("p"), ...)
                 panel.loess(x, y, col=4)
                 panel.grid(lty=3, col="lightgrey")   })
xyplot(y ~ temp | as.factor(time), 
       data = dialectric, col=1,
       panel = function(x,y,...) 
               { panel.xyplot(x, y, type=c("p"), ...)
                 panel.loess(x, y, col=4) 
                 panel.grid(lty=3, col="lightgrey")  })

dialectric$log.time <- log(dialectric$time)
mod.glm <- glm(y ~ log.time + temp, data = dialectric, family=Gamma(inverse))
summary(mod.glm)
plot(mod.glm)

mod1 <- fwdglm(y ~ log.time + temp, data = dialectric, 
               family=Gamma(inverse), nsamp=500)
plot(mod1, 1, th.Res = 1)
plot(mod1, 8)

#------------------------------------------------#
# Script for computing Goodness Of Link test 
# at specified lambda values

lambda <- c("inverse", "1/mu^2", "log", "sqrt", "identity", "mu^2")
link.test <- matrix(NA, length(lambda), 1)
links <- NULL
for (j in 1:length(lambda))
    { if (lambda[j]=="mu^2")
         { q <- quasi(link="identity", var="mu^2")
           q[c("linkfun", "linkinv", "mu.eta", "valideta")] <- power(2)
           q$link <- lambda[j]
         }  
      else   
         { q <- quasi(link=lambda[j], var="mu^2") }
      links[[j]] <- q                                    
    }

for (j in 1:length(lambda))    
    { g <- glm(y ~ log.time + temp, data = dialectric, family = links[[j]], x=T, y=T)
      disp <- sum(residuals(g, "pearson")^2)/g$df.residual
      link.test[j,1] <- scglm(g$x, g$y, family = g$family, weights = g$weights, 
                              beta = g$coefficients, phi = disp)
    }
rownames(link.test) <- lambda; colnames(link.test) <- "Link test"
print(link.test)

mod2 <- fwdglm(y ~ log.time + temp, data = dialectric, family=links[[6]], nsamp=500)
plot(mod2, 8)
plot(mod2, 1, th.Res = 1)

mod3.glm <- glm(y ~ log.time*temp + I(log.time^2) + I(temp^2), 
                data = dialectric, family=Gamma(inverse), data = dialectric)

mod3 <- fwdglm(y ~ log.time*temp + I(log.time^2) + I(temp^2), 
               data = dialectric, family=Gamma(inverse),  nsamp=500)
plot(mod3,1)

D1 <- rep(0,length(dialectric$y)); D1[c(109:112,125:128)] <- 1 
D2 <- rep(0,length(dialectric$y)); D2[c(93:96)] <- 1 

col <- pch <- rep(1,length(D1))
col[D1==1] <- pch[D1==1] <- 2
col[D2==1] <- pch[D2==1] <- 3
xyplot(y ~ log.time | as.factor(temp), data = dialectric,
       groups=as.factor(temp), col=col, pch=pch,
       panel = function(x,y,subscripts,col,pch,...)
               { panel.xyplot(x, y, type=c("p"), 
                              col=col[subscripts], 
                              pch=pch[subscripts], ...)
                 panel.loess(x, y, col=1)  
                 panel.grid(lty=3, col="lightgrey") })

mod4.glm <- glm(y ~ log.time*temp + I(log.time^2) + I(temp^2) + D1 + D2, 
                data = dialectric, family=Gamma(inverse),)
summary(mod4.glm)
               
mod4 <- fwdglm(y ~ log.time*temp + I(log.time^2) + I(temp^2) + D1 + D2, 
               data = dialectric, family=Gamma(inverse), nsamp=500)
plot(mod4, 1, th.Res=0.5)

D3 <- rep(0,length(dialectric$y)); D3[c(61:64,77:80)] <- 1 

col <- pch <- rep(1,length(dialectric$y))
col[D1==1] <- pch[D1==1] <- 2
col[D2==1] <- pch[D2==1] <- 3
col[D3==1] <- pch[D3==1] <- 4
xyplot(y ~ log.time | as.factor(temp), data = dialectric,
       groups=as.factor(temp), col=col, pch=pch,
       panel = function(x,y,subscripts,col,pch,...)
               { panel.xyplot(x, y, type=c("p"), 
                              col=col[subscripts], 
                              pch=pch[subscripts], ...)
                 panel.loess(x, y, col=1)  
                 panel.grid(lty=3, col="lightgrey") })

mod5.glm <- glm(y ~ log.time*temp + I(log.time^2) + I(temp^2) + D1 + D2 + D3,
                data = dialectric, family=Gamma(inverse))
summary(mod5.glm)
anova(mod5.glm, test="F")


#------------------------------------------------#
# Script for running the forward-search and 
# compute the goodness of link test at specified 
# lambda values

lambda <- c("inverse", "1/mu^2", "log", "sqrt", "identity", "mu^2")
mod5 <- NULL
# create a list with link functions
links <- NULL
for (j in 1:length(lambda))
    { if (lambda[j]=="mu^2")
         { q <- quasi(link="identity", var="mu^2")
           q[c("linkfun", "linkinv", "mu.eta", "valideta")] <- power(2)
           q$link <- lambda[j]
         }  
      else   
         { q <- quasi(link=lambda[j], var="mu^2") }
      links[[j]] <- q                                    
    }
# run the forward search for each link function separately 
# and save the results in the list mod5
for (j in 1:length(lambda))    
    { g <- fwdglm(y ~ log.time*temp + I(log.time^2) + I(temp^2) + D1 + D2 + D3,
                  data = dialectric, family = links[[j]], nsamp=500)
      mod5[[j]] <- g
      names(mod5)[[j]] <- paste("Link =", lambda[[j]])
    }
# plot the result
for (j in 1:length(mod5))
    { 
      m <- sapply(mod5[[j]]$included, length)[-1]
      if (j==1)
         plot(mod5[[j]], 8, ylim=c(-5,5), xlim=c(5.5, 140))
      else
         lines(m, mod5[[j]]$ScoreTest, lty=j, col=j)
      text(max(m)+1, mod5[[j]]$ScoreTest[diff(range(m))+1], lambda[j], pos=4)
    }
#------------------------------------------------#

plot(mod5$"Link = sqrt", 5)

mod6 <- fwdglm(y ~ log.time*temp + I(log.time^2) + D1 + D2 + D3,
               family = links[[4]], data = dialectric, nsamp=1000)
summary(mod6)
plot(mod6, 1, th.Res=0.3, ylim=c(-0.8, 0.6))
plot(mod6, 2)
plot(mod6, 5)
plot(mod6, 6)
plot(mod6, 7)
plot(mod6, 8, ylim=c(-3,3))
plot(mod6, 3)


#-----------------------------------------------------------------------------#

data(derailme)
derailme$Type <- as.factor(derailme$Type)

library(lattice)
trellis.device(bg="white")

xyplot(log(y/TrainKm) ~ Year | Type, 
       data = derailme, col=1, layout=c(3,1), 
       strip = function(...) 
               strip.default(..., strip.names=c(T,T), style=3),
       panel = function(x,y,...) 
               { panel.xyplot(x, y, type=c("p"), ...)
                 panel.loess(x, y, col=4)
                 panel.grid(lty=3, col="lightgrey")  })

mod.glm <- glm(y ~ Year + Type + offset(log(TrainKm)), 
               data = derailme, family=poisson(log))
summary(mod.glm)

j <- c(13, 23, 63)  # exclude the three largest obs
mod.glm <- glm(y ~ Year + Type + offset(log(TrainKm)), 
               data = derailme, family=poisson(log), subset=-j)
summary(mod.glm)

plot(fitted(mod.glm), residuals(mod.glm), type="n")
text(fitted(mod.glm), residuals(mod.glm), derailme$y[-j], 
     col=as.numeric(derailme$Type)[-j])
lim <- tapply(fitted(mod.glm), derailme$Type[-j], max)
abline(v=c(1.568, 2.5), lty=3)
text(1.2, -1.7, adj=0, "Non pass.")
text(1.8, -1.7, adj=0, "Post mark 1")
text(2.6, -1.7, adj=0, "Mark 1")

mod1 <- fwdglm(y ~ Year + factor(Type) + offset(log(TrainKm)),
               data = derailme, family=poisson(log), nsamp = 1000)
summary(mod1)
plot(mod1, 6)

mod2 <- fwdglm(y ~ factor(Type) + offset(log(TrainKm)),
               data = derailme, family=poisson(log), nsamp = 1000)
summary(mod2)
plot(mod2, 1)
plot(mod2, 8)
plot(mod2, 7)

dev <- mod2$Likelihood[,1]
p <- ncol(mod2$Coefficients)
m <- sapply(mod2$included, length)
plot(m, dev, type="l", xlab="Subset Size", ylab="Deviance")
axis(1, at=seq(0,max(m), by=5)); grid()
lines(m, qchisq(0.95, m-3), lty=2, col=2)
lines(m, qchisq(0.99, m-3), lty=3, col=2)

#-----------------------------------------------------------------------------#

data(cellular)
cellular$tnf <- cellular$TNF
cellular$ifn <- cellular$IFN
cellular$TNF <- as.factor(cellular$TNF)
cellular$IFN <- as.factor(cellular$IFN)

mod.glm <- glm(y ~ TNF + IFN, family=poisson(log), data = cellular, x=T, y=T)
summary(mod.glm)
plot(mod.glm)

mod1 <- fwdglm(y ~ TNF + IFN, family=poisson(log), data = cellular, nsamp=200)
summary(mod1)

plot(mod1, 1, th.Res=1.5)
plot(mod1, 7)
plot(mod1, 10)
plot(mod1, 8)

cellular$TNF <- as.numeric(cellular$TNF)
cellular$IFN <- as.factor(cellular$IFN)

library(lattice)
trellis.device(bg="white")

xyplot(y ~ log(tnf+1) | factor(log(ifn+1)), 
       data = cellular, col=1,  
       strip = function(...) 
               strip.default(..., strip.names=c(T,T), style=3),
       panel = function(x,y,...) 
               { panel.xyplot(x, y, type=c("b"), ...)
                 panel.grid(lty=3, col="lightgrey")  })

xyplot(y ~ log(ifn+1) | factor(log(tnf+1)), 
       data = cellular, col=1,  
       strip = function(...) 
               strip.default(..., strip.names=c(T,T), style=3),
       panel = function(x,y,...) 
               { panel.xyplot(x, y, type=c("b"), ...)
                 panel.grid(lty=3, col="lightgrey")  })


mod2 <- fwdglm(y[-16] ~ (log(tnf+1)[-16])*(log(tnf+1)[-16]),
               family=poisson(log), data = cellular, nsamp=200)
plot(mod2, 7)
plot(mod2, 6)
plot(mod2, 8)

#-----------------------------------------------------------------------------#

data(bliss)
bliss$log.dose <- log(bliss$Dose) 

mod1 <- fwdglm(y ~ log.dose, weights=Total, data = bliss,
               family=binomial(logit), nsamp="all")
summary(mod1)

mod2 <- fwdglm(y ~ log.dose, weights=Total, data = bliss,
               family=binomial(probit), nsamp="all")
summary(mod2)

mod3 <- fwdglm(y ~ log.dose, weights=Total, data = bliss,
               family=binomial(cloglog), nsamp="all")
summary(mod3)

oldpar <- par(); par(mfrow=c(3,1), mar=c(4,4,3,2))
plot(mod1, 1, squared=T, th.Res=2.5); title("Logit link")
plot(mod2, 1, squared=T, th.Res=1.5); title("Probit link")
plot(mod3, 1, squared=T, th.Res=1); title("Cloglog link")
par(oldpar)

plot(mod1, 8, ylim=c(-3,3))
m <- sapply(mod1$included, length)
lines(m[-1], mod2$ScoreTest, lty=2)
lines(m[-1], mod3$ScoreTest, lty=3)
legend(3,-2, c("Logit link", "Probit link", "Cloglog link"), lty=1:3)

x0 <- seq(3, 5, length=100)
# logit link
inc <- mod1$included
m2 <- glm(y ~ log.dose, weights=Total, data = bliss,
          subset=inc[[1]], family=binomial(logit))
m8 <- glm(y ~ log.dose, weights=Total, data = bliss,
          subset=inc[[7]], family=binomial(logit))
fit.logit.m2 <- predict(m2, data.frame(log.dose=x0), type="response")
fit.logit.m8 <- predict(m8, data.frame(log.dose=x0), type="response")
# cloglog link
inc <- mod2$included
m2 <- glm(y ~ log.dose, weights=Total, data = bliss, 
          subset=inc[[1]], family=binomial(cloglog))
m8 <- glm(y ~ log.dose, weights=Total, data = bliss, 
          subset=inc[[7]], family=binomial(cloglog))
fit.cloglog.m2 <- predict(m2, data.frame(log.dose=x0), type="response")
fit.cloglog.m8 <- predict(m8, data.frame(log.dose=x0), type="response")
#
oldpar <- par(); par(mfrow=c(2,2), mar=c(4,4,3,2))
#
plot(y~log.dose, data=bliss, xlim=range(x0), ylim=c(0,1))
lines(x0, fit.logit.m2)
title("Logit link (m=2)")
text(bliss$log.dose, bliss$y, rownames(bliss), cex=0.7, pos=4, col=3)
#
plot(y~log.dose, data=bliss, xlim=range(x0), ylim=c(0,1))
lines(x0, fit.logit.m8)
title("Logit link (m=8)")
text(bliss$log.dose, bliss$y, rownames(bliss), cex=0.7, pos=4, col=3)
#
plot(y~log.dose, data=bliss, xlim=range(x0), ylim=c(0,1))
lines(x0, fit.cloglog.m2)
title("Cloglog link (m=2)")
text(bliss$log.dose, bliss$y, rownames(bliss), cex=0.7, pos=4, col=3)
#
plot(y~log.dose, data=bliss, xlim=range(x0), ylim=c(0,1))
lines(x0, fit.cloglog.m8)
title("Cloglog link (m=8)")
text(bliss$log.dose, bliss$y, rownames(bliss), cex=0.7, pos=4, col=3)
par(oldpar)

#-----------------------------------------------------------------------------#

data(mice)
mice$y <- mice$conv/mice$total
mice$log.dose <- log(mice$dose)

mark <- mice$prep+1
plot(y ~ log.dose, data = mice, ylim=c(0,1), col=mark, pch=mark)
text(mice$log.dose, mice$y, rownames(mice), cex=0.7, col=mark, pos=4, offset=0.3)

mod1.glm <- glm(y ~ log.dose + prep, weights=total, data = mice, 
               family=binomial(logit)) 
summary(mod1.glm)
mod2.glm <- glm(y ~ log.dose + prep, weights=total, data = mice, 
                family=binomial(cloglog))
summary(mod2.glm)
mod3.glm <- glm((1-y) ~ log.dose + prep, weights=total, data = mice, 
                family=binomial(cloglog))
summary(mod3.glm)

mod1 <- fwdglm(y ~ log.dose + prep, weights=total, data = mice,
               family=binomial(logit), nsamp="all") 
plot(mod1, 1, th.Res=1, ylim=c(-3,3))

mod2 <- fwdglm(y ~ log.dose + prep, weights=total, data = mice,
               family=binomial(cloglog), nsamp="all") 
plot(mod2, 1, th.Res=0.8, ylim=c(-3,3))

mod3 <- fwdglm((1-y) ~ log.dose + prep, weights=total, data = mice, 
               family=binomial(cloglog), nsamp="all") 
plot(mod3, 1, th.Res=1, ylim=c(-3,3))

plot(mod1, 8, ylim=c(-3,3))
m <- sapply(mod1$included, length)
lines(m[-1], mod2$ScoreTest, lty=2)
lines(m[-1], mod3$ScoreTest, lty=3)
legend(4,-1, c("Logit link", "Cloglog link", "Loglog link"), lty=1:3)

m <- sapply(mod1$included, length)
plot(0,0, xlim=range(m), ylim=c(0,14), xlab="Subset Size", ylab="Deviance", type="n")
lines(m, mod1$Likelihood[,1], lty=1)
lines(m, mod2$Likelihood[,1], lty=2)
lines(m, mod3$Likelihood[,1], lty=3)
legend(4,14, c("Logit link", "Cloglog link", "Loglog link"), lty=1:3)

#-----------------------------------------------------------------------------#

data(rainfall)
rainfall$z <- scale(rainfall$Rain)
rainfall$z2 <- rainfall$z^2
rainfall$z3 <- rainfall$z^3
rainfall$y <- rainfall$Cases/rainfall$Total

mod0.glm <- glm(y ~ 1, weights=Total, data = rainfall, 
                family=binomial(logit))
mod1.glm <- glm(y ~ z + I(z^2) + I(z^3), weights=Total, data = rainfall, 
                family=binomial(logit)) 
summary(mod1.glm)
1-pchisq(62.635, 30)
anova(mod0.glm, mod1.glm, test="Chisq")

# Very sensible to carefully choosen initial subset. 
mod1 <- fwdglm(y ~ z + z2 + z3, weights=Total, data = rainfall, 
               family=binomial(logit), nsamp="all") 
# which takes very long, and the best starting subset obtained is
inibsb <- c(1, 2, 3, 28)    
mod1 <- fwdglm(y ~ z + z2 + z3, weights=Total, data = rainfall, 
               family=binomial(logit), bsb=inibsb) 

plot(mod1, 1, th.Res=2)
plot(mod1, 8)
plot(mod1, 7)
plot(mod1, 9)
plot(mod1, 6)

inc <- mod1$included
inc.m34 <- inc$"m=34"
inc.m30 <- inc$"m=30"
inc.m29 <- inc$"m=29"

plot(y ~ z, data = rainfall, ylim=c(0,1))
j <- setdiff(inc.m34, inc.m30)
text(rainfall$z[j], rainfall$y[j], j, col=2, cex=0.7, offset=0.3, pos=4)
j <- setdiff(inc.m30, inc.m29)
text(rainfall$z[j], rainfall$y[j], j, col=3, cex=0.7, offset=0.3, pos=4)
z0 <- seq(min(rainfall$z), max(rainfall$z), length=100)
lines(z0, predict(mod1.glm, data.frame(z=z0), type="response"))
lines(z0, predict(update(mod1.glm, subset=inc.m30), data.frame(z=z0), type="response"), lty=2, col=2)
lines(z0, predict(update(mod1.glm, subset=inc.m29), data.frame(z=z0), type="response"), lty=3, col=3)
legend(0,1,c("m=34", "m=30", "m=29"), col=1:3, lty=1:3)

# Function for computing deletion residuals
"delres.glm" <- function(obj)
{ 
  r <- residuals(obj, type="pearson")
  disp <- summary(obj)$dispersion
  lev <- influence.measures(obj)$infmat[,8]
  r/(disp*sqrt(1-lev))
}

q <- qqnorm(delres.glm(mod1.glm))
text(q$x[abs(q$y)>2], q$y[abs(q$y)>2], which(abs(q$y)>2), 
     cex=0.7, offset=0.3, pos=4)

#-----------------------------------------------------------------------------#

data(vaso)
vaso$log.vol <- log(vaso$volume)
vaso$log.rate <- log(vaso$rate)

col <- pch <- rep(1, length(vaso$y))
col[vaso$y==1] <- 2; pch[vaso$y==1] <- 4
j <- c(4,18,24,17,32)
plot(log.rate ~ log.vol, data = vaso, pch=pch, col=col)
text(vaso$log.vol[j], vaso$log.rate[j], j, pos=4, offset=0.3, cex=0.8)

# logistic model on the original data
mod.glm <- glm(y ~ log.vol + log.rate, data = vaso, family=binomial(logit), x=T)
summary(mod.glm)

r <- residuals(mod.glm)
fit.linpred <- predict(mod.glm)
fit.prob <- predict(mod.glm, type="response")

plot(fit.prob, r, pch=pch, col=col)
text(fit.prob[j], r[j], j, pos=4, offset=0.3, cex=0.8)

plot(y ~ fit.linpred, data = vaso, pch=pch, col=col, ylim=c(0,1.05))
lines(sort(fit.linpred), 1/(1+exp(-sort(fit.linpred))))
text(fit.linpred[j], vaso$y[j], j, pos=3, offset=0.3, cex=0.8)

# logistic model on the modified data
out <- c(4,18,24)
mod1.glm <- glm(y ~ log.vol + log.rate, data = vaso, 
                family=binomial(logit), subset=-out, x=T)
summary(mod1.glm)

fit.linpred <- predict(mod1.glm, data.frame(mod.glm$x[,-1]))
plot(y ~ fit.linpred, data = vaso, pch=pch, col=col, ylim=c(0,1.05))
lines(sort(fit.linpred), 1/(1+exp(-sort(fit.linpred))))
text(fit.linpred[out], vaso$y[out], out, pos=3, offset=0.3, cex=0.8)

# quadratic logistic model on the original data
mod2.glm <- glm(y ~ log.vol*log.rate+I(log.vol^2)+I(log.rate^2), data = vaso,
                family=binomial(logit))
anova(mod2.glm, test="Chisq")

# balanced search
mod1 <- fwdglm(y ~ log.vol + log.rate, data = vaso, 
               family=binomial(logit), nsamp=1000)

plot(mod1, 6, ylim=c(-4,4), plot.pf=TRUE)
plot(mod1, 1, ylim=c(-2.5,8), xlim=c(36,39.2), th.Res=8, plot.pf=TRUE)

# unbalanced search
mod2 <- fwdglm(y ~ log.vol + log.rate, data = vaso, 
              family=binomial(logit), balanced=F,  nsamp=1000)


#-----------------------------------------------------------------------------#

data(chapman)

mark <- chapman$y+1
pairs(chapman[,c(1,4,6)], col=mark, pch=mark)

mod1 <- fwdglm(y ~ age + chol + weight, data = chapman, 
               family=binomial(logit), epsilon=10^-8, nsamp=1000)
plot(mod1, 6, plot.pf=TRUE)
plot(mod1, 7, plot.pf=TRUE)

