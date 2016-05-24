### lmrob()  with "real data"

library(robustbase)

set.seed(0)
data(salinity)
summary(m0.sali  <- lmrob(Y ~ . , data = salinity))
(A1 <- anova(m0.sali, Y ~ X1 + X3))
## -> X2 is not needed
(m1.sali  <- lmrob(Y ~ X1 + X3, data = salinity))
(A2 <- anova(m0.sali, m1.sali)) # the same as before
stopifnot(all.equal(A1[2,"Pr(>chisq)"],
		    A2[2,"Pr(>chisq)"], tolerance=1e-14))
anova(m0.sali, m1.sali, test = "Deviance")
## whereas 'X3' is highly significant:
m2 <- update(m0.sali, ~ . -X3)
(A3 <- anova(m0.sali, m2))
(A4 <- anova(m0.sali, m2, test = "Deviance"))
cX3 <- c(Estimate = -0.627327396, `Std. Error` = 0.15844971,
         `t value` = -3.9591577, `Pr(>|t|)` = 0.000584156)
stopifnot(all.equal(cX3, coef(summary(m0.sali))["X3",], tolerance = 1e-6))


##  example(lmrob)
set.seed(7)
data(coleman)
summary( m1 <- lmrob(Y ~ ., data=coleman) )
stopifnot(c(3,18) == which(m1$w < 0.2))

data(starsCYG)
(RlmST <- lmrob(log.light ~ log.Te, data = starsCYG))
summary(RlmST)
stopifnot(c(11,20,30,34) == which(RlmST$w < 0.01))

set.seed(47)
data(hbk)
m.hbk <- lmrob(Y ~ ., data = hbk)
summary(m.hbk)
stopifnot(1:10 == which(m.hbk$w < 0.01))

data(heart)
summary(mhrt <- lmrob(clength ~ ., data = heart))
stopifnot(8  == which(mhrt$w < 0.15),
          11 == which(0.61 < mhrt$w & mhrt$w < 0.62),
          c(1:7,9:10,12) == which(mhrt$w > 0.90))

data(stackloss)
mSL <- lmrob(stack.loss ~ ., data = stackloss)
summary(mSL)


cat('Time elapsed: ', proc.time(),'\n') # "stats"
