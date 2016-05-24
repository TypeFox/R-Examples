### tests methods argument of lmrob.control

library(robustbase)

data(stackloss)

## S
set.seed(0)
summary(m0 <- lmrob(stack.loss ~ ., data = stackloss, method = "S",
                    compute.outlier.stats = "S"))
set.seed(0)
m0a <- lmrob.S(m0$x, stack.loss, lmrob.control())

all.equal(m0[c('coefficients', 'scale', 'rweights')],
          m0a[c('coefficients', 'scale', 'rweights')])

## MM
set.seed(0)
summary(m1 <- lmrob(stack.loss ~ ., data = stackloss, method = "MM",
                    compute.outlier.stats = "S"))

set.seed(0)
m2 <- update(m1, method = "SM")

all.equal(m1[c('coefficients', 'scale', 'cov')],
          m2[c('coefficients', 'scale', 'cov')])

set.seed(0)
m3 <- update(m0, method = "SM", cov = '.vcov.w')

## SMD
set.seed(0)
summary(m4 <- lmrob(stack.loss ~ ., data = stackloss, method = "SMD", psi = 'bisquare',
                    compute.outlier.stats = "S"))
summary(m4a <- lmrob..D..fit(m3))

## rearrange m4a and update call
m4a <- m4a[names(m4)]
class(m4a) <- class(m4)
m4a$call <- m4$call

all.equal(m4, m4a)

## SMDM
set.seed(0)
summary(m5 <- lmrob(stack.loss ~ ., data = stackloss, method = "SMDM", psi = 'bisquare',
                    compute.outlier.stats = "S"))
summary(m5a <- lmrob..M..fit(obj=m4))

## rearrange m5a
m5a <- m5a[names(m5)]
class(m5a) <- class(m5)

all.equal(m5, m5a)

## Fast S large n strategy (sped up)
model <- model.frame(LNOx ~ . ,data = NOxEmissions)
control <- lmrob.control(fast.s.large.n = 10, n.group = 341, groups = 2)
set.seed(0)
try(ret <- lmrob.S(model.matrix(model, NOxEmissions)[1:682,], NOxEmissions$LNOx[1:682], control))
## do what the error says
control <- lmrob.control(fast.s.large.n = Inf)
try(ret <- lmrob.S(model.matrix(model, NOxEmissions)[1:682,], NOxEmissions$LNOx[1:682], control))
## this still fails, but this error is to be expected since only a part
## of the design matrix is given
