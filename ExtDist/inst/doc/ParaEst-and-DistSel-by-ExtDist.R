## ---- message=FALSE------------------------------------------------------
library(ExtDist)

## ----Data, comment=""----------------------------------------------------
set.seed(1234)
head(X <- rWeibull(50, shape = 2, scale = 3))

## ----Weibull est, comment=""---------------------------------------------
est.par <- eWeibull(X)

## ----class, comment=""---------------------------------------------------
class(est.par)

## ----dpqr egs, results='hide'--------------------------------------------
dWeibull(seq(0,2,0.4), params = est.par)
pWeibull(seq(0,2,0.4), params = est.par)
qWeibull(seq(0,1,0.2), params = est.par)
rWeibull(10, params = est.par)

## ---- results='hide'-----------------------------------------------------
dWeibull(seq(0,2,0.4), shape = est.par$shape, scale = est.par$scale)
pWeibull(seq(0,2,0.4), shape = est.par$shape, scale = est.par$scale)
qWeibull(seq(0,1,0.2), shape = est.par$shape, scale = est.par$scale)
rWeibull(10, shape = est.par$shape, scale = est.par$scale)

## ----selection criterion, comment=""-------------------------------------
logLik(est.par) # log likelihood
AIC(est.par) # Akaike information criterion
AICc(est.par) # corrected Akaike information criterion
BIC(est.par) # Bayes' Information Criterion. 
MDL(est.par) # minimum description length 
vcov(est.par) # variance-covariance matrix of the parameters of the fitted distribution

## ----Example, comment=""-------------------------------------------------
Ozone <- airquality$Ozone
Ozone <- Ozone[!is.na(Ozone)] # Removing the NA's from Ozone data
summary(Ozone)
best <- bestDist(Ozone, candDist=c("Gamma", "Weibull", "Normal", "Exp"), criterion = "logLik");best

## ----DistSelCriteria, comment=""-----------------------------------------
DistSelCriteria(Ozone, candDist = c("Gamma", "Weibull", "Normal", "Exp"),
                         criteria = c("logLik","AIC","AICc", "BIC"))

## ----CompareDist, comment=""---------------------------------------------
compareDist(Ozone, attributes(best)$best.dist.par, eNormal(Ozone))

## ----ploteDist, error=FALSE----------------------------------------------
plot(attributes(best)$best.dist.par)

## ----Chunk10-------------------------------------------------------------
Y <- c(0.1703, 0.4307, 0.6085, 0.0503, 0.4625, 0.479, 0.2695, 0.2744, 0.2713, 0.2177, 
       0.2865, 0.2009, 0.2359, 0.3877, 0.5799, 0.3537, 0.2805, 0.2144, 0.2261, 0.4016)
w <- c(0.85, 1.11, 0.88, 1.34, 1.01, 0.96, 0.86, 1.34, 0.87, 1.34, 0.84, 0.84, 0.83, 1.09, 
       0.95, 0.77, 0.96, 1.24, 0.78, 1.12)

## ----Chunk11, comment=""-------------------------------------------------
eBeta(Y,w)

bestDist(Y, w, candDist = c("Beta_ab","Laplace","Normal"), criterion = "AIC")

DistSelCriteria(Y, w, candDist = c("Beta_ab","Laplace","Normal"),
                         criteria = c("logLik","AIC","AICc", "BIC"))

