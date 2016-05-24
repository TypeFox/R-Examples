## Example of observation-specific boundaries
## Estimate the willingness to pay for the Kakadu National Park
## Data given in intervals -- 'lower' for lower bound and 'upper' for upper bound.
## Note that dichotomous-coice answers are already coded to 'lower' and 'upper'
set.seed(1)
options(digits=4)
library(intReg)
data(Kakadu, package="Ecdat")
## Estimate in log form, change 999 to Inf
Kakadu <- Kakadu[sample(nrow(Kakadu), 300),]
                           # Speed up the tests
lb <- log(Kakadu$lower)
ub <- Kakadu$upper
ub[ub > 998] <- Inf
ub <- log(ub)
y <- cbind(lb, ub)
m <- intReg(y ~ sex + log(income) + age + schooling + 
              recparks + jobs + lowrisk + wildlife + future + aboriginal + finben +
              mineparks + moreparks + gov +
              envcon + vparks + tvenv + major, data=Kakadu)
## You may want to compare the results to Werner (1999),
## Journal of Business and Economics Statistics 17(4), pp 479-486

## Test coef, stdEr, summary with and without boundaries
print(coef(m))
print(coef(m, boundaries=TRUE))
print(stdEr(m))
print(stdEr(m, boundaries=TRUE))
print(summary(m))
print(summary(m, boundaries=TRUE))

## test model.matrix
mm <- model.matrix(m)
print(mm[i <- sample(nrow(mm), 10),])

##
## Example of common intervals for all the observations
##
library(Ecdat)
data(Bwages)
## calculate an ordinary Mincer-style wage regression.  
## Note: gross hourly wage rate in EUR
intervals <- c(0, 5, 10, 15, 25, Inf)
salary <- cut(Bwages$wage, intervals)
int <- intReg(salary ~ factor(educ) + poly(exper, 2), data=Bwages, boundaries=log(intervals))
## Note: use logs for the intervals in Euros.  We do not have to
## transform salaris to log form as this does not change the intervals.
## Ignore any warnings
cat("Interval regression:\n")
print(summary(int))

## Test model.frame
mf <- model.frame(int)
print(mf[i <- sample(nrow(mf), 10),])
print(model.response(mf)[i])

##
## Small data, large number of intervals (by Thierry Kalisa)
##
a <- c(0.002300, 0.020000, 0.000150, 0.000005, 0.002300, 0.000045, 0.000150,
       0.000110, 0.000110, 0.000005, 0.010000, 0.000490, 0.000110, 0.000005,
       0.000600, 0.000380, 0.000600, 0.005275, 0.005275, 0.000045, 0.000075,
       0.000600, 0.000600, 0.005275, 0.000075, 0.001650, 0.001100, 0.000005,
       0.000025, 0.005275, 0.000150, 0.005275, 0.000005, 0.000110, 0.000270,
       0.000600, 0.000600, 0.000380, 0.000110, 0.000380, 0.000270, 0.000490,
       0.000045, 0.000110, 0.000110, 0.000150, 0.000005, 0.000110, 0.000045,
       0.005275, 0.000600, 0.000200, 0.003475, 0.005275, 0.000005, 0.000600,
       0.000200, 0.000075, 0.000600, 0.000600, 0.000075, 0.000230, 0.000490,
       0.005275, 0.000230, 0.000110, 0.000490, 0.000045, 0.000075, 0.001650,
       0.000600, 0.000490, 0.000005, 0.003475, 0.001650, 0.000150, 0.000380,
       0.017500, 0.003475, 0.000270, 0.000230, 0.005275, 0.000045, 0.000045,
       0.000075, 0.003475, 0.000150, 0.002300, 0.001650, 0.001100, 0.000005,
       0.000075, 0.000025, 0.000025, 0.000150, 0.001100)
b <- c(0.003475, 0.040000, 0.005275, 0.040000, 0.015000, 0.001100, 0.000380,
       0.003475, 0.003475, 0.040000, 0.020000, 0.007075, 0.000490, 0.003475,
       0.007075, 0.005275, 0.012500, 0.012500, 0.010000, 0.000270, 0.000200,
       0.002300, 0.010000, 0.010000, 0.001650, 0.003475, 0.005275, 0.003475,
       0.003475, 0.010000, 0.000600, 0.020000, 0.000045, 0.001650, 0.010000,
       0.005275, 0.020000, 0.001650, 0.005275, 0.003475, 0.003475, 0.007075,
       0.002300, 0.010000, 0.000270, 0.000270, 0.003475, 0.000600, 0.000270,
       0.007075, 0.003475, 0.010000, 0.010000, 0.012500, 0.000045, 0.010000,
       0.003475, 0.010000, 0.012500, 0.003475, 0.000380, 0.003475, 0.005275,
       0.008650, 0.000600, 0.002300, 0.003475, 0.005275, 0.003475, 0.003475,
       0.003475, 0.002300, 0.000025, 0.017500, 0.005275, 0.003475, 0.001650,
       0.020000, 0.040000, 0.001650, 0.003475, 0.008650, 0.000200, 0.000110,
       0.000490, 0.040000, 0.000600, 0.020000, 0.005275, 0.008650, 0.000490,
       0.005275, 0.000230, 0.000200, 0.000270, 0.005275)
c <-c(3, 4, 3, 3, 3, 1, 2, 1, 3, 4, 2, 2, 1, 2, 1, 2, 2, 1, 3, 2, 2, 3, 1, 2, 1, 2, 3, 2, 4, 3, 4, 2,
      4, 2, 1, 2, 4, 3, 2, 3, 2, 2, 3, 4, 2, 1, 3, 3, 1, 1, 2, 1, 2, 2, 1, 3, 1, 1, 2, 3, 2, 2, 3, 1,
      3, 2, 2, 1, 2, 2, 2, 2, 1, 3, 2, 3, 2, 1, 1, 2, 2, 1, 1, 2, 3,
      1, 2, 3, 2, 2, 1, 1, 4, 1, 3, 3)
ab <- cbind(a,b)
mNorm <- intReg(ab~c)
print(summary(mNorm))

## Test the same with cloglog disturbances
mCloglog <- intReg(ab~c, method="cloglog")
