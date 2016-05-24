## Example of ordinary least squares

x <- rep(1,15)
y <- x + rnorm(15)

mod1 <- lm(y~1)
mod2 <- regress(y~1)

## Test if coefficients are the same
if(abs(coef(mod1) - mod2$beta) > 1e-10) stop("Error with Ordinary Least Squares with constant")

## Test estimate of sigma-squared
if(abs(anova(mod1)[,3] - mod2$sigma) > 1e-10) stop("Error with Ordinary Least Squares with constant (variance)")


## Examples from lm

## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.  
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2,10,20, labels=c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)
lm.D90 <- lm(weight ~ group - 1) # omitting intercept

rg.D9 <- regress(weight ~ group)
rg.D90 <- regress(weight ~ group - 1)

if( sum((coef(lm.D9) - rg.D9$beta)^2 + (coef(lm.D90) - rg.D90$beta)^2 > 1e-10)) stop("Error with Ordinary Least Squares")

if(abs(anova(lm.D9)[2,3] - rg.D9$sigma) > 1e-10) stop("Error with Ordinary Least Squares (variance)")
if(abs(anova(lm.D90)[2,3] - rg.D90$sigma) > 1e-10) stop("Error with Ordinary Least Squares (variance)")
