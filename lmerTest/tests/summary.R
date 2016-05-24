require(lmerTest)

load(system.file("testdata", "ProblemSet.RData", package="lmerTest"))
m <- lmer(ExpScore ~ ChLevel+(1|Subject), data=ProblemSet)
summary(m)

## bug from Cyrus
test=data.frame(x = rnorm(100),f=factor(rep(c(1:10),10)))
test.mod=lmer(x~1+(1|f),data=test)
summary(test.mod)


require(lmerTest)

carrots$Income <- as.factor(carrots$Income)
carrots$Age <- as.factor(carrots$Age)
carrots$Frequency <- as.factor(carrots$Frequency)
carrots$Gender <- as.factor(carrots$Gender)

l <- list(Frequency="contr.SAS", Income="contr.SAS")
m.carrots <- lmer(Preference ~ sens2*Frequency*Income
                  +(1+sens2|Consumer), data=carrots, contrasts=l)
an.m <- anova(m.carrots)

TOL <- 1e-4 # for the check

#with 4 decimals should agree with SAS output
#numbers before decimals should agree with SAS output
stopifnot(
  all.equal(an.m[,"Pr(>F)"], c(2e-5, 0.15512,  0.06939, 0.08223, 0.52459, 0.03119, 0.48344), tol = TOL), 
  all.equal(round(an.m$DenDF), c(83, 83, 83, 83, 83, 83, 83))
  , TRUE)



sm <- summary(m.carrots)
stopifnot(
  all.equal(sm$coefficients[,"Pr(>|t|)"], c(1e-10, 0.005061, 0.6865554, 0.342613, 0.129157, 
                                            0.088231, 0.846000, 0.354472, 0.526318, 0.020646, 0.010188, 
                                            0.031242, 0.055356, 0.694689, 0.099382, 0.28547, 
                                            0.977774, 0.855653, 0.427737, 0.321086, 0.417465 , 0.204385, 0.784437,
                                            0.681434, 0.106180, 0.149122, 0.390870, 0.273686), tol=TOL,
            check.attributes = FALSE), TRUE)

sm.kr <- summary(m.carrots, ddf = "Kenward-Roger")

## coefficients for Sat and KR agree in this example
all.equal(sm$coefficients[,"Pr(>|t|)"], sm$coefficients[,"Pr(>|t|)"])


## checking lsmeans and difflsmeans
## compare with SAS output
m <- lmer(Informed.liking ~ Product*Information*Gender 
          + (1|Product:Consumer) + (1|Consumer) , data=ham)


lsm <- lsmeans(m, test.effs = "Product")

TOL <- 1e-7
stopifnot(
  all.equal( lsm$lsmeans.table[, "Estimate"], c(5.8084, 5.1012, 6.0909, 5.9256),             
             tol=TOL,
             check.attributes = FALSE),            
  all.equal(lsm$lsmeans.table[, "t-value"], c(24.93, 21.89, 26.14, 25.43), tol=TOL,
            check.attributes = FALSE),
  all.equal(lsm$lsmeans.table[, "Lower CI"], c(5.3499, 4.6428, 5.6324, 5.4672), tol=TOL,
            check.attributes = FALSE),
  all.equal(lsm$lsmeans.table[, "Upper CI"], c(6.2668, 5.5597, 6.5493, 6.3840), tol=TOL,
            check.attributes = FALSE),
  TRUE)



