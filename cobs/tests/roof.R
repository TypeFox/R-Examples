suppressMessages(library(cobs))

data(USArmyRoofs)
attach(USArmyRoofs)#-> "age" and "fci"

if(!dev.interactive(orNone=TRUE)) pdf("roof.pdf", width=10)

## Compute the quadratic median smoothing B-spline with SIC
## chosen lambda
a50 <- cobs(age,fci,constraint = "decrease",lambda = -1,nknots = 10,
            degree = 2,pointwise = rbind(c(0,0,100)),
            trace = 2)# trace > 1 : more tracing
summary(a50)

## Generate Figure 2a (p.22 in online version)
## Plot $pp.lambda againt $pp.sic :
plot(a50$pp.lambda, a50$pp.sic, type = "b", log = "x",
     xlab = "Lambda", ylab = "SIC") ##<< no lambda in [20, 180] -- bug !?
lam0 <- a50$pp.lambda[6] ## should be the 2nd smallest one
abline(v = c(a50$lambda, lam0),lty = 2)

## For "testing" (signif() to stay comparable):
signif(cbind(a50$pp.lambda, a50$pp.sic), 6)

## Compute the quadratic median smoothing B-spline with lambda
## at the the 2nd smallest SIC
a50.1 <- cobs(age,fci,constraint = "decrease",lambda = lam0, nknots = 10,
              degree = 2,pointwise = rbind(c(0,0,100)))
summary(a50.1)
## Compute the 25th percentile smoothing B-spline
a25 <- cobs(age,fci,constraint = "decrease",lambda = -1, nknots = 10,
            tau = .25, pointwise = rbind(c(0,0,100)))
summary(a25)

## Again plot $pp.sic against $pp.lambda
plot(a25$pp.lambda,a25$pp.sic,type = "l",log = "x")
abline(v = a25$lambda,lty = 2)

## Compute the 75th percentile smoothing B-spline
a75 <- cobs(age, fci, constraint = "decrease",lambda = -1, nknots = 10,
            tau = .75, pointwise = rbind(c(0, 0, 100)))
a75
zapsmall(a75$coef)

## We rerun cobs with ... as suggested by warning of scobs (? right ??)
## a75 <- cobs(age, fci, constraint = "decrease", lambda = -1, nknots = 10,
##             tau = .75,pointwise = rbind(c(0, 0, 100)))
## ## still changing tau , ok ?
## a75 <- cobs(age, fci, constraint = "decrease", lambda = -1, nknots = 10,
##             tau = .75,pointwise = rbind(c(0, 0, 100)))
## summary(a75)

## Again we plot $pp.sic against $pp.lambda
plot(a75$pp.lam,a75$pp.sic,type = "l",log = "x")
## It seems like the linear fit is really what the data wants

(pa50   <- predict(a50,  interval = "both"))
(pa50.1 <- predict(a50.1,interval = "both"))
(pa25   <- predict(a25,  interval = "both"))
(pa75   <- predict(a75,  interval = "both"))
## Generate Figure 2b (p.22, online version)
plot(age,fci,xlim = c(0,15),ylim = c(0,100),xlab = "AGE",ylab = "FCI")
lines(pa50)
lines(pa50.1,lty = 2)
lines(pa25, col = 2)
lines(pa75, col = 3)

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
