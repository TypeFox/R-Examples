### Demonstrate 'treatReg' usage

DGP <- function(N=1000, sigma=1, rho=0.8,
                alpha0=-1, alpha1=1, alpha2=1,
                beta0=0, beta1=1, beta2=1) {
   ## Generate random data
   library(mvtnorm)
   Sigma <- matrix(c(1, rho*sigma, rho*sigma, sigma^2), 2, 2)
   uv <- rmvnorm(N, mean=c(0,0), sigma=Sigma)
   u <- uv[,1]
   v <- uv[,2]
   x <- rnorm(N)
   z <- rnorm(N)
   ySX <- alpha0 + alpha1*x + alpha2*z + u
   yS <- ySX > 0
   yO <- beta0 + beta1*x + beta2*yS + v
   data.frame(yO, yS, x, z, ySX, u, v)
}

library(sampleSelection)
## Create random data:
## yS = 1(-1 + x + z + u > 0
## yO = x + yS + v
## x, z are standard normals
## (u,v) ~ N((0,0), matrix(c(1, 0.8, 0.8, 1), 2, 2))
## -- the correlation between disturbance terms 0.8
## -- the true treatment effect is 1.0
dat <- DGP(2000)
## 1) Estimate the treatment effect by OLS:
ols <- lm(yO ~ x + yS, data=dat)
print(summary(ols))
## The effect (yS) is substantially overestimated (ySTRUE close to 2)
readline("<press Enter>")
## Now estimate the same model with treatReg:
tr <- treatReg(yS~x+z, yO~x+yS, data=dat)
print(summary(tr))
## -- now the effect (yS) is close to 1
## next: a real example

readline("<press Enter>")
data(Treatment, package="Ecdat")
## 'Treatment' data in library 'Ecdat':
## 'treat'   treatment indicator (logical)
## 'age'     in years
## 'educ'    education in years
## 'u74', 'u75'  unemployment in 1974, 1975 (logical)
## 'ethn'    race: black, hispanic, other
## 're78'    real income in 1978
## First estimate it using unemployment 'u74', 'u75' as exclusion
## restriction
er <- treatReg(treat~poly(age,2) + educ + u74 + u75 + ethn,
               log(re78)~treat + poly(age,2) + educ + ethn,
               data=Treatment)
print(summary(er))
## The treatment effect estimate 'treatTRUE' is -0.96, i.e. the
## treatment substantially lower the earnings
## Now estimate it withouth the exclusion restriction
noer <- treatReg(treat~poly(age,2) + educ + u74 + u75 + ethn,
                 log(re78)~treat + poly(age,2) + educ + u74 + u75 + ethn,
                 data=Treatment)
print(summary(noer))
## Now the estimate is -0.51. The treatment still seems to
## lower the earnings.
## The few observable characteristics are probably not sufficient
## to correct for selection effects.  Model selection and
## valid exclusion restrictions are crucial in this type of models

