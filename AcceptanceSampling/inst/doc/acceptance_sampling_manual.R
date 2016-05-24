### R code from vignette source 'acceptance_sampling_manual.Rnw'

###################################################
### code chunk number 1: acceptance_sampling_manual.Rnw:605-606
###################################################
library(AcceptanceSampling)


###################################################
### code chunk number 2: acceptance_sampling_manual.Rnw:622-624
###################################################
x <- OC2c(10, 3)
x


###################################################
### code chunk number 3: acceptance_sampling_manual.Rnw:629-631
###################################################
plot(x)
grid(lty="solid")


###################################################
### code chunk number 4: acceptance_sampling_manual.Rnw:653-660
###################################################
xb <- OC2c(5, 1, type="b")                     ## Binomial
xh <- OC2c(5, 1, type="h", N=50, pd=(0:50)/50) ## Hypergeometric
xp <- OC2c(5, 1, type="p")                     ## Poisson
plot(xb, type="l", xlim=c(0, 0.2), ylim=c(0.6, 1))
grid(lty="solid")
points(xh@pd, xh@paccept, col="green")
lines(xp@pd, xp@paccept, col="red")


###################################################
### code chunk number 5: acceptance_sampling_manual.Rnw:694-699
###################################################
x.mean <- seq(248, 255, 0.05)
x.pd <- pnorm(250, mean=x.mean, sd=1.5)
x.plan <- OC2c(10, 1, pd=x.pd)
plot(x.mean, x.plan, xlab="Mean weight")
grid(lty="solid")


###################################################
### code chunk number 6: acceptance_sampling_manual.Rnw:715-717
###################################################
x <- OC2c(10,3, pd=seq(0,0.1,0.01))
summary(x, full=TRUE)


###################################################
### code chunk number 7: acceptance_sampling_manual.Rnw:733-734
###################################################
assess(OC2c(20,0), PRP=c(0.05, 0.95), CRP=c(0.15, 0.075))


###################################################
### code chunk number 8: acceptance_sampling_manual.Rnw:749-750
###################################################
find.plan(PRP=c(0.05, 0.95), CRP=c(0.15, 0.075), type="binom")


###################################################
### code chunk number 9: acceptance_sampling_manual.Rnw:786-790
###################################################
x <- OC2c(n=c(8,8), c=c(0,1), r=c(2,2))
x
plot(x)
grid(lty="solid")


###################################################
### code chunk number 10: acceptance_sampling_manual.Rnw:825-831
###################################################
x.mean <- seq(248, 255, 0.05)
x.pd <- pnorm(250, mean=x.mean, sd=1.5)
find.plan(PRP=c(0.05, 0.95), CRP=c(0.15, 0.075), type="normal", s.type="known")
x.plan <- OCvar(n=26, k=1.322271, pd=x.pd)
plot(x.mean, x.plan, xlab="Mean weight")
grid(lty="solid")


###################################################
### code chunk number 11: acceptance_sampling_manual.Rnw:857-858
###################################################
find.plan(PRP=c(0.05, 0.95), CRP=c(0.15, 0.075), type="normal", s.type="unknown")


###################################################
### code chunk number 12: acceptance_sampling_manual.Rnw:866-867
###################################################
find.plan(PRP=c(0.05, 0.95), CRP=c(0.15, 0.075), type="normal", s.type="known")


###################################################
### code chunk number 13: acceptance_sampling_manual.Rnw:883-890
###################################################
xb <- OC2c(n=80,c=7)
xn1 <- OCvar(n=49, k=1.326538, s.type="unknown")
xn2 <- OCvar(n=26, k=1.322271)
plot(xb, type="l", xlim=c(0,0.3))
grid(lty="solid")
lines(xn1@pd, xn1@paccept, col="green")
lines(xn2@pd, xn2@paccept, col="red")


