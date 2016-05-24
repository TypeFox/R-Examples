library(cobs)
data(USArmyRoofs)
attach(USArmyRoofs)#-> "age" and "fci"

if(!interactive()) postscript("roof.ps", horizontal = TRUE)

## Compute the quadratic median smoothing B-spline with SIC
## chosen lambda
a50 <- cobs(age,fci,constraint = "decrease",lambda = -1,nknots = 10,
            degree = 2,pointwise = rbind(c(0,0,100)),
            lstart = 1e3, trace = 2)# trace > 1 : more tracing

## Generate Figure 2a (p.22 in online version)
## Plot $pp.lambda againt $sic :
plot(a50$pp.lambda, a50$sic, type = "b", log = "x",
     xlab = "Lambda", ylab = "SIC") ##<< no lambda in [20, 180] -- bug !?
lam0 <- a50$pp.lambda[6] ## should be the 2nd smallest one
abline(v = c(a50$lambda, lam0),lty = 2)

## For "testing" (signif() to stay comparable):
signif(cbind(a50$pp.lambda, a50$sic), 6)

## Compute the quadratic median smoothing B-spline with lambda
## at the the 2nd smallest SIC
a50.1 <- cobs(age,fci,constraint = "decrease",lambda = lam0, nknots = 10,
              degree = 2,pointwise = rbind(c(0,0,100)), lstart = 1e3)
## Compute the 25th percentile smoothing B-spline
a25 <- cobs(age,fci,constraint = "decrease",lambda = -1, nknots = 10,
            tau = .25, pointwise = rbind(c(0,0,100)))
## Again plot $sic against $pp.lambda
plot(a25$pp.lambda,a25$sic,type = "l",log = "x")
abline(v = a25$lambda,lty = 2)

## Compute the 75th percentile smoothing B-spline
a75 <- cobs(age, fci, constraint = "decrease",lambda = -1, nknots = 10,
            tau = .75,pointwise = rbind(c(0, 0, 100)))
## We rerun cobs with a larger lstart at 10^8 as suggested by warning of cobs
a75 <- cobs(age, fci, constraint = "decrease",lambda = -1, nknots = 10,
            tau = .75,lstart = 10^8,pointwise = rbind(c(0, 0, 100)))
## still lstart...
a75 <- cobs(age, fci, constraint = "decrease",lambda = -1, nknots = 10,
            tau = .75,lstart = 10^10,pointwise = rbind(c(0, 0, 100)))

## Again we plot $sic against $pp.lambda
plot(a75$pp.lam,a75$sic,type = "l",log = "x")
## It seems like the linear fit is really what the data wants

## Generate Figure 2b (p.22, online version)
plot(age,fci,xlim = c(0,15),ylim = c(0,100),xlab = "AGE",ylab = "FCI")
pa50 <- predict(a50,interval = "none")
pa50.1 <- predict(a50.1,interval = "none")
pa25 <- predict(a25,interval = "none")
pa75 <- predict(a75,interval = "none")
lines(pa50[,1],pa50[,2])
lines(pa50.1[,1],pa50.1[,2],lty = 2)
lines(pa25[,1],pa25[,2], col = 2)
lines(pa75[,1],pa75[,2], col = 3)
