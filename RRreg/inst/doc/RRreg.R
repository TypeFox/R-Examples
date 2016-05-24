## ------------------------------------------------------------------------
library(RRreg)
data.W <- RRgen(n=1000, pi.true=.3, model="Warner", p=.2, complyRates = c(1, 1), sysBias = c(0,0))
head(data.W)

## ------------------------------------------------------------------------
warner <- RRuni(response=response, data=data.W, model="Warner", p=.2, MLest=T)
summary(warner)

## ------------------------------------------------------------------------
data.W$cov[data.W$true==1] <- rnorm(sum(data.W$true == 1),1)
data.W$cov[data.W$true==0] <- rnorm(sum(data.W$true == 0))

## ------------------------------------------------------------------------
RRcor(x=data.W$response, y=data.W$cov, models=c("Warner", "direct"), p.list=list(.2), bs.n=0, bs.type = c("se.n", "se.p", "pval"), nCPU=1)

## ------------------------------------------------------------------------
cor(x=data.W$true, y=data.W$cov)

## ------------------------------------------------------------------------
log1 <- RRlog(formula=response~cov, data=data.W, model="Warner", p=.2, LR.test=TRUE, fit.n = 1)
summary(log1)

## ------------------------------------------------------------------------
glm(formula=true~cov, data=data.W, family=binomial(link = "logit"))

## ------------------------------------------------------------------------
# make random cluster:
data.W$cluster <- c("a", "b")
mixmod <- RRmixed(response ~ cov + (1 | cluster), data=data.W, model = "Warner", p = .2)
mixmod

## ------------------------------------------------------------------------
lin1 <- RRlin(formula=cov~response, data=data.W, models="Warner", p.list=.2, bs.n=0, fit.n=1)
summary(lin1)

## ------------------------------------------------------------------------
data.W$pred <- rnorm(1000)
lin2 <- RRlin(formula=cov~response + pred, data=data.W, models="Warner", p.list=list(.2), bs.n=0, fit.n=1)
summary(lin2)

## ------------------------------------------------------------------------
lm(cov~true + pred, data=data.W)

## ------------------------------------------------------------------------
data.W2 <- RRgen(n=1000, pi.true=.45, model="Warner", p=.35)
data.W$cov <-2*data.W$true - 3*data.W2$true+ rnorm(1000,1,5)
data.W$response2 <- data.W2$response
lin3 <- RRlin(cov ~ response + response2, data=data.W, models=c("Warner", "Warner"), p.list=list(.2, .35), fit.n=1)
summary(lin3)

## ------------------------------------------------------------------------
# generate data with different prevalence rates for RR and DQ
RR <- RRgen(n=500, pi=.4, model="Warner", p=.35)
DQ <- rbinom(500, 1, .2)
response <- c(RR$response, DQ)
# dummy variable for RR vs. DQ
group <- rep(1:0, each=500)
# include interaction of question format and age (here, the rescaled z-variable z.age)
z.age <- c(rnorm(500, mean=RR$true*2,sd=3), rnorm(500, mean=0,sd=3))
# fit full model (i.e., test difference in prevalence estimates and interaction)
fit <- RRlog(response ~ group * z.age, model="Warner", p=.35, group=group, LR.test=TRUE, fit.n = 1) 
summary(fit)
# get prevalence estimate for RR and DQ
logit1 <- coef(fit) %*% c(1, 1, mean(z.age), mean(z.age))
cat(exp(logit1) / (1+ exp(logit1)))
logit2 <- coef(fit) %*% c(1, 0, mean(z.age), 0)
cat(exp(logit2) / (1+ exp(logit2)))

## ------------------------------------------------------------------------
getPW(model = "FR", p = c(.1, .1))

## ------------------------------------------------------------------------
### generate 2 two-group RR variables and a continuous covariate
RR1 <- RRgen(1000, pi=.4, model="SLD", p=c(.2,.8), complyRates=c(.8,1))
RR2 <- RRgen(1000, pi=.6, model="SLD", p=c(.3 ,.7), complyRates=c(.9,1))
cov <- RR1$true + RR2$true + rnorm(1000, 0, 3)

### logistic regression (dependent RR variable)
logmod <- RRlog(RR1$response ~ cov, model="SLD", p=c(.2,.8), group=RR1$group, fit.n = 1)
summary(logmod)

### group matrix for RRlin and RRcor: use multiple group vectors in one matrix
group <- cbind(RR1$group, RR2$group)

# bivariate correlation between 2 two-group variables and a continuous covariate
# note that only the most informative subsets of the data are used (see ?RRcor)
responses <- cbind(RR1$response, RR2$response, cov)
colnames(responses) <- c("RR1", "RR2", "cov")
RRcor( responses, models=c("SLD","SLD","direct"), p.list=list(c(.2,.8), c(.3, .7)), group=group)

# linear model with 2 RR predictors
linmod <- RRlin(cov ~ RR1$response+RR2$response, models=c("SLD","SLD"), p.list=list(c(.2,.8), c(.3, .7)), group=group, fit.n=1)
summary(linmod)

