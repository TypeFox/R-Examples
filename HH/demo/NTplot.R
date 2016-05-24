## 22 distinct cases of NTplot (right-sided shown).  All these also work with the additional argument shiny=TRUE.

NTplot(mean1=2, xbar=1.5, xlim=c(-3, 5), n=4                                            )
NTplot(mean1=2, xbar=1.5, xlim=c(-3, 5), n=4, number.vars=2                             )
NTplot(mean1=2, xbar=1.5, xlim=c(-3, 5), n=4,                distribution.name="t", df=3)
NTplot(mean1=2, xbar=1.5, xlim=c(-3, 5), n=4, number.vars=2, distribution.name="t", df=3)

NTplot(mean1=2, xbar=1.5, xlim=c(-3, 5), n=4                                            , type="confidence")
NTplot(mean1=2, xbar=1.5, xlim=c(-3, 5), n=4, number.vars=2                             , type="confidence")
NTplot(mean1=2, xbar=1.5, xlim=c(-3, 5), n=4,                distribution.name="t", df=3, type="confidence")
NTplot(mean1=2, xbar=1.5, xlim=c(-3, 5), n=4, number.vars=2, distribution.name="t", df=3, type="confidence")


x1 <- c(5,3,2,7,6,0,4,3)
x2 <- c(7,3,5,9,4,3,7,4)
t1 <- t.test(x1, mu=2, alt="greater")
t1
t2 <- t.test(x2, x1, mu=-1, alt="greater")
t2
t2p <- t.test(x2, x1, mu=-1, alt="greater", paired=TRUE)
t2p
t2ps <- t.test(x2-x1, mu=-1, alt="greater")  ## verify paired
t2ps

NTplot(t1, mean1=5)
NTplot(t1, type="confidence")
NTplot(t2, mean1=2)
NTplot(t2, type="confidence")
NTplot(t2p, mean1=1.6)
NTplot(t2p, type="confidence")

pt1 <-  power.t.test(power = .90, delta = 1.2, sd=2, type="one.sample", alternative = "one.sided")
pt1
pt2 <-  power.t.test(power = .90, delta = 1.2, sd=2,                    alternative = "one.sided")
pt2
pt2p <- power.t.test(power = .90, delta = 1.2, sd=2, type="paired",     alternative = "one.sided")
pt2p

NTplot(pt1,  xbar=1)
NTplot(pt1,  xbar=1, type="confidence")
NTplot(pt2,  xbar=1)
NTplot(pt2,  xbar=1, type="confidence")
NTplot(pt2p, xbar=1)
NTplot(pt2p, xbar=1, type="confidence")

NTplot(p0=.4, p.hat=.65, p1=.7, distribution.name="binomial", n=15)
NTplot(p.hat=.65, distribution.name="binomial", n=15, type="confidence")


## display options
NTplot(t1, mean1=5)
NTplot(t1, mean1=5, power=TRUE)
NTplot(t1, mean1=5, power=TRUE, beta=TRUE)
NTplot(t1, mean1=5,             beta=TRUE)
print(NTplot(t1, mean1=5), tablesOnPlot=FALSE)
print(NTplot(t1, mean1=5, power=TRUE, beta=TRUE), tablesOnPlot=FALSE)




## Power Calculations for Two-Sample Test for Proportions
PPT <- power.prop.test(n = 50, p1 = .50, p2 = .75, alternative = "one.sided")
PPT
try(NTplot(PPT)) ## not yet working


## Test of Equal or Given Proportions
## From ?prop.test
     ## Data from Fleiss (1981), p. 139.
     ## H0: The null hypothesis is that the four populations from which
     ##     the patients were drawn have the same true proportion of smokers.
     ## A:  The alternative is that this proportion is different in at
     ##     least one of the populations.

     smokers  <- c( 83, 90, 129, 70 )
     patients <- c( 86, 93, 136, 82 )
     prop.test(smokers, patients)  ## four groups
PT2 <- prop.test(smokers[3:4], patients[3:4], alternative = "greater")  ## two groups
PT2
try(NTplot(PT2)) ## not yet working
PT1 <- prop.test(smokers[4], patients[4], p=.75, alternative = "greater")  ## one group
PT1
try(NTplot(PT1)) ## not yet working
