### R code from vignette source 'simulation.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: simulation.Rnw:53-55
###################################################
library(actuar)
options(width = 62, digits = 4)


###################################################
### code chunk number 2: simulation.Rnw:89-90 (eval = FALSE)
###################################################
## rpois(n, rgamma(n, 3, rgamma(n, 2, 2)))


###################################################
### code chunk number 3: simulation.Rnw:266-267
###################################################
set.seed(3)


###################################################
### code chunk number 4: simulation.Rnw:269-280
###################################################
nodes <- list(cohort = 2,
              contract = c(4, 3),
              year = c(4, 4, 4, 4, 5, 5, 5))
mf <- expression(cohort = rexp(2),
                 contract = rgamma(cohort, 1),
                 year = rpois(weights * contract))
ms <- expression(cohort = rnorm(2, sqrt(0.1)),
                 contract = rnorm(cohort, 1),
                 year = rlnorm(contract, 1))
wijt <- runif(31, 0.5, 2.5)
pf <- simul(nodes = nodes, model.freq = mf, model.sev = ms, weights = wijt)


###################################################
### code chunk number 5: simulation.Rnw:287-290
###################################################
class(pf)
pf$data
pf$classification


###################################################
### code chunk number 6: simulation.Rnw:302-303
###################################################
pf


###################################################
### code chunk number 7: simulation.Rnw:311-313
###################################################
aggregate(pf)
aggregate(pf, by = c("cohort", "year"), FUN = mean)


###################################################
### code chunk number 8: simulation.Rnw:320-322
###################################################
frequency(pf)
frequency(pf, by = "cohort")


###################################################
### code chunk number 9: simulation.Rnw:338-340
###################################################
severity(pf)
severity(pf, splitcol = 1)


###################################################
### code chunk number 10: simulation.Rnw:345-346
###################################################
weights(pf)


###################################################
### code chunk number 11: simulation.Rnw:351-352
###################################################
aggregate(pf, classif = FALSE) / weights(pf, classif = FALSE)


###################################################
### code chunk number 12: simulation.Rnw:380-382
###################################################
set.seed(123)
options(width = 55)


###################################################
### code chunk number 13: simulation.Rnw:384-386
###################################################
wit <- rgamma(15, rep(runif(3, 0, 100), each = 5),
              rep(runif(3, 0, 100), each = 5))


###################################################
### code chunk number 14: simulation.Rnw:394-398
###################################################
frequency(simul(list(entity = 3, year = 5),
      expression(entity = rgamma(rgamma(1, 5, 5), rgamma(1, 25, 1)),
          year = rpois(weights * entity)),
      weights = wit))


