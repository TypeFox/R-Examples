### R code from vignette source 'Devore7.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(width=75, show.signif.stars = FALSE)
library(Devore7)
library(lattice)


###################################################
### code chunk number 2: medexample
###################################################
conc = c(7.6, 8.3, 9.3, 9.4, 9.4, 9.7, 10.4, 11.5, 11.9, 15.2, 16.2, 20.4)
str(conc)
median(conc)


###################################################
### code chunk number 3: libraryD6
###################################################
library(Devore7)


###################################################
### code chunk number 4: dataxmp0113
###################################################
data(xmp01.13)
str(xmp01.13)


###################################################
### code chunk number 5: withxmp
###################################################
with(xmp01.13, median(concentration))


###################################################
### code chunk number 6: Devore7.Rnw:166-167 (eval = FALSE)
###################################################
## library(Devore7)


###################################################
### code chunk number 7: Devore7.Rnw:174-176 (eval = FALSE)
###################################################
## data(xmp01.13)
## str(xmp01.13)


###################################################
### code chunk number 8: xmp01.01
###################################################
with(xmp01.01, stem(temp))


###################################################
### code chunk number 9: Devore7.Rnw:251-252
###################################################
with(xmp01.01, hist(temp))


###################################################
### code chunk number 10: xmp01.05
###################################################
with(xmp01.05, stem(bingePct))
with(xmp01.05, stem(bingePct, scale = 0.5))


###################################################
### code chunk number 11: xmp01.06
###################################################
with(xmp01.06, stem(yardage))


###################################################
### code chunk number 12: xmp0109
###################################################
with(xmp01.09, hist(consump))


###################################################
### code chunk number 13: xmp0111
###################################################
with(xmp01.10, hist(strength, breaks = c(2,4,6,8,12,20,30)))


###################################################
### code chunk number 14: xmp0112
###################################################
with(xmp01.12, mean(crackLength))
with(xmp01.12, sum(crackLength))
with(xmp01.13, median(concentration))
with(xmp01.12, summary(crackLength))
with(xmp01.13, summary(concentration))


###################################################
### code chunk number 15: xmp0114
###################################################
with(xmp01.14, summary(copper))
with(xmp01.14, mean(copper, trim = 0.1))


###################################################
### code chunk number 16: xmp0115
###################################################
with(xmp01.15, var(Strength))
with(xmp01.15, sd(Strength))


###################################################
### code chunk number 17: setup1
###################################################
opar = par(mar=c(3.5,4.1,0.1,0.1))


###################################################
### code chunk number 18: xmp0117
###################################################
with(xmp01.17, boxplot(depth, horizontal = TRUE))


###################################################
### code chunk number 19: setup2
###################################################
par(opar)


###################################################
### code chunk number 20: setup1
###################################################
opar = par(mar=c(3.5,4.1,0.1,0.1))


###################################################
### code chunk number 21: xmp0118
###################################################
with(xmp01.18, boxplot(C1, horizontal = TRUE))


###################################################
### code chunk number 22: setup2
###################################################
par(opar)


###################################################
### code chunk number 23: setup1
###################################################
opar = par(mar=c(3.5,4.1,0.1,0.1))


###################################################
### code chunk number 24: xmp0119
###################################################
with(ex01.15, boxplot(Score ~ Type, horizontal = TRUE, las = 1))


###################################################
### code chunk number 25: setup2
###################################################
par(opar)


###################################################
### code chunk number 26: Devore7.Rnw:476-477
###################################################
with(xmp01.09, hist(consump, las = 1))


###################################################
### code chunk number 27: Devore7.Rnw:486-487
###################################################
with(xmp01.01, stem(temp, scale=2))


###################################################
### code chunk number 28: Devore7.Rnw:490-491
###################################################
with(xmp01.01, hist(temp, breaks=c(25,35,45,55,65,75,85)))


###################################################
### code chunk number 29: xmp0223
###################################################
choose(15,3)*choose(10,3)/choose(25,6)


###################################################
### code chunk number 30: xmp0330
###################################################
dbinom(3, size = 6, prob = 0.5)
dbinom(3:6, size = 6, prob = 0.5)
sum(dbinom(3:6, size = 6, prob = 0.5))
1 - pbinom(2, size = 6, prob = 0.5)
pbinom(2, size = 6, prob = 0.5, lower = FALSE)
pbinom(1, 6, 0.5)


###################################################
### code chunk number 31: xmp0331
###################################################
pbinom(8,15,0.2)
dbinom(8,15,0.2)
1-pbinom(7,15,0.2)
sum(dbinom(4:7,15,0.2))


###################################################
### code chunk number 32: xmp0335
###################################################
dhyper(2,12,8,5)


###################################################
### code chunk number 33: xmp0336
###################################################
dhyper(2,5,20,10)
phyper(2,5,20,10)


###################################################
### code chunk number 34: xmp0338
###################################################
dnbinom(10, 5, 0.2)
pnbinom(10, 5, 0.2)


###################################################
### code chunk number 35: xmp0339
###################################################
dpois(5, lambda = 4.5)
ppois(5, 4.5)


###################################################
### code chunk number 36: xmp0340
###################################################
dbinom(1, 400, 0.005)
dpois(1,2)
pbinom(3, 400, 0.005)
ppois(3, 2)


###################################################
### code chunk number 37: xmp0413a
###################################################
pnorm(1.25)
pnorm(-1.25)


###################################################
### code chunk number 38: xmp0413b
###################################################
1-pnorm(1.25)
pnorm(1.25, lower=FALSE)


###################################################
### code chunk number 39: xmp0413d
###################################################
pnorm(c(-0.38,1.25))
diff(pnorm(c(-0.38,1.25)))


###################################################
### code chunk number 40: xmp0414
###################################################
qnorm(0.99)


###################################################
### code chunk number 41: xmp0415
###################################################
qnorm(0.05, lower=FALSE)


###################################################
### code chunk number 42: xmp0416
###################################################
diff(pnorm(c(1.0, 1.75), mean = 1.25, sd = 0.46))


###################################################
### code chunk number 43: xmp0418
###################################################
qnorm(0.995, mean=64, sd=0.78)


###################################################
### code chunk number 44: xmp0420
###################################################
pnorm(10.5, mean = 12.5, sd = sqrt(12.5*0.75))
pbinom(10, size = 50, prob = 0.25)
diff(pnorm(c(4.5,15.5), mean = 12.5, sd = sqrt(12.5*0.75)))
diff(pbinom(c(4,15), size = 50, prob = 0.25))


###################################################
### code chunk number 45: xmp0423
###################################################
pgamma(c(3,5), shape=2)
diff(pgamma(c(3,5), shape=2))
pgamma(4, shape = 2, lower = FALSE)


###################################################
### code chunk number 46: xmp0424
###################################################
diff(pgamma(c(60,120), shape = 8, scale = 15))
pgamma(30, shape = 8, scale = 15, lower = FALSE)


###################################################
### code chunk number 47: xmp0421
###################################################
pexp(10, rate = 0.2)
diff(pexp(c(5,10), rate = 0.2))
pexp(2, rate = 0.5, lower = FALSE)


###################################################
### code chunk number 48: xmp0425
###################################################
pweibull(10, shape = 2, scale = 10)
qweibull(0.95, shape = 2, scale = 10)


###################################################
### code chunk number 49: xmp0427
###################################################
diff(plnorm(c(1,2), meanlog = 0.375, sdlog = 0.25))
qlnorm(0.99, meanlog = 0.375, sdlog = 0.25)


###################################################
### code chunk number 50: xmp0427
###################################################
pbeta((3-2)/(5-2), shape1 = 2, shape2 = 3)


###################################################
### code chunk number 51: xmp0429
###################################################
with(xmp04.29, qqnorm(meas.err))
with(xmp04.29, qqline(meas.err))


###################################################
### code chunk number 52: xmp0430
###################################################
with(xmp04.30, qqnorm(Voltage))
with(xmp04.30, qqline(Voltage))


###################################################
### code chunk number 53: xmp0431
###################################################
with(xmp04.31, plot(log(-log(1-seq(0.05, 0.95, 0.1))), log(lifetime)))


###################################################
### code chunk number 54: xmp0522aa
###################################################
options(digits=5)


###################################################
### code chunk number 55: xmp0522a
###################################################
curve(dweibull(x, shape = 2, scale = 5), 0, 15, las = 1)


###################################################
### code chunk number 56: xmp0522b
###################################################
rweibull(10, shape = 2, scale = 5)


###################################################
### code chunk number 57: xmp0522c
###################################################
samp = rweibull(10, shape = 2, scale = 5)
print(samp)
mean(samp)
median(samp)
sd(samp)


###################################################
### code chunk number 58: xmp0522d
###################################################
means = medians = sds = numeric(6)
for (i in 1:6) {
    samp = rweibull(10, shape = 2, scale = 5)
    print(samp)
    means[i] = mean(samp)
    medians[i] = median(samp)
    sds[i] = sd(samp)
}
means
medians
sds


###################################################
### code chunk number 59: xmp0522zz
###################################################
options(digits=7)


###################################################
### code chunk number 60: xmp0523a
###################################################
par(mfrow = c(2,2))
samp5 = matrix(rnorm(500 * 5, mean = 8.25, sd = 0.75), ncol = 500)
hist(colMeans(samp5), main='Samples of size 5')
samp10 = matrix(rnorm(500 * 10, mean = 8.25, sd = 0.75), ncol = 500)
hist(colMeans(samp10), main='Samples of size 10')
samp20 = matrix(rnorm(500 * 20, mean = 8.25, sd = 0.75), ncol = 500)
hist(colMeans(samp20), main='Samples of size 20')
samp30 = matrix(rnorm(500 * 30, mean = 8.25, sd = 0.75), ncol = 500)
hist(colMeans(samp30), main='Samples of size 30')


###################################################
### code chunk number 61: xmp0523z
###################################################
par(mfrow=c(1,1))


###################################################
### code chunk number 62: xmp0523a
###################################################
curve(dlnorm(x, meanlog=3, sdlog=0.4), from = 0, to = 75, las = 1)


###################################################
### code chunk number 63: xmp0523b
###################################################
par(mfrow=c(2,2))
samp5=matrix(rlnorm(500 * 5, 2, 0.4), ncol = 500)
hist(colMeans(samp5), main = 'Means of samples of size 5')
samp10=matrix(rlnorm(500 * 10, 2, 0.4), ncol = 500)
hist(colMeans(samp10), main = 'Means of samples of size 10')
samp20=matrix(rlnorm(500 * 20, 2, 0.4), ncol = 500)
hist(colMeans(samp20), main = 'Means of samples of size 20')
samp30=matrix(rlnorm(500 * 30, 2, 0.4), ncol = 500)
hist(colMeans(samp30), main = 'Means of samples of size 30')


###################################################
### code chunk number 64: xmp0523c
###################################################
qqnorm(colMeans(samp30))
qqline(colMeans(samp30))


###################################################
### code chunk number 65: xmp0602
###################################################
with(xmp06.02, mean(Voltage))
with(xmp06.02, median(Voltage))
with(xmp06.02, mean(range(Voltage)))
with(xmp06.02, mean(Voltage, trim=0.1))


###################################################
### code chunk number 66: xmp0603a
###################################################
with(xmp06.03, var(Strength))
with(xmp06.03, sd(Strength))


###################################################
### code chunk number 67: xmp0603b
###################################################
with(xmp06.03, sum((Strength-mean(Strength))^2)/length(Strength))


###################################################
### code chunk number 68: xmp0613a
###################################################
xbar = with(xmp06.13, mean(Survival))
xsqb = with(xmp06.13, mean(Survival^2))
xbar
xsqb
xbar^2/(xsqb-xbar^2)
(xsqb-xbar^2)/xbar


###################################################
### code chunk number 69: xmp0612b
###################################################
library(MASS)
with(xmp06.13, fitdistr(Survival, dgamma, list(shape=10.577,scale=10.726)))


###################################################
### code chunk number 70: xmp0702
###################################################
80.0 + c(-1, 1) * 1.96 * 2.0 / sqrt(31)


###################################################
### code chunk number 71: xmp0706
###################################################
with(xmp07.06, mean(Voltage)+c(-1,1)*1.96*sd(Voltage)/sqrt(length(Voltage)))


###################################################
### code chunk number 72: xmp0706b
###################################################
with(xmp07.06, t.test(Voltage))


###################################################
### code chunk number 73: xmp0706c
###################################################
with(xmp07.06, qqnorm(Voltage))


###################################################
### code chunk number 74: xmp0708
###################################################
prop.test(16,48)
binom.test(16,48)


###################################################
### code chunk number 75: xmp0711a
###################################################
with(xmp07.11, t.test(Elasticity))
with(xmp07.11, qqnorm(Elasticity))


###################################################
### code chunk number 76: xmp0715
###################################################
with(xmp07.15,qqnorm(voltage))
with(xmp07.15, 16*var(voltage)/qchisq(c(0.975,0.025), df = 16))


###################################################
### code chunk number 77: xmp0808
###################################################
with(xmp08.08, t.test(DCP, mu = 30, alt = "less"))


###################################################
### code chunk number 78: xmp0808b
###################################################
with(xmp08.08, qqnorm(DCP))


###################################################
### code chunk number 79: xmp0808c
###################################################
with(xmp08.08, qqnorm(log(DCP)))


###################################################
### code chunk number 80: xmp08.08d
###################################################
with(xmp08.08, t.test(log(DCP), mu=log(30), alt="less"))


###################################################
### code chunk number 81: xmp0809
###################################################
with(xmp08.09, t.test(MAWL, mu = 25, alt = "greater"))


###################################################
### code chunk number 82: xmp0810
###################################################
power.t.test(n=10,delta=0.1,sd=0.1,type="one.sample",alt="one.sided")
power.t.test(delta=0.1,sd=0.1,power=0.95,type="one.sample",alt="one.sided")


###################################################
### code chunk number 83: xmp0811a
###################################################
prop.test(1276, 4115, p = 0.3, alt = "greater")


###################################################
### code chunk number 84: xmp0811b
###################################################
prop.test(1276, 4115, p = 0.3, alt = "greater", correct=FALSE)
sqrt(1.993)


###################################################
### code chunk number 85: xmp0813
###################################################
binom.test(x=14,n=20,p=0.9,alt="less")


###################################################
### code chunk number 86: strbw
###################################################
print(bwplot(type ~ strength, xmp09.07, xlab = "Tensile strength (psi)"))


###################################################
### code chunk number 87: strqq
###################################################
print(qqmath(~ strength|type, xmp09.07,
      xlab = "Standard normal quantiles",
      ylab = "Tensile strength (psi)", 
      type = c("g","p"), aspect = 1))


###################################################
### code chunk number 88: xmp0907
###################################################
str(xmp09.07)


###################################################
### code chunk number 89: xmp0907bwprt (eval = FALSE)
###################################################
## bwplot(type ~ strength, xmp09.07)


###################################################
### code chunk number 90: xmp0907qqprt (eval = FALSE)
###################################################
## qqmath(~strength|type, xmp09.07)


###################################################
### code chunk number 91: xmp0907t
###################################################
t.test(strength ~ type, xmp09.07, alt = "greater")


###################################################
### code chunk number 92: zincdiffqq
###################################################
print(qqmath( ~ I(bottom-surface), xmp09.08, type = c("g","p"),
      xlab = "Standard normal quantiles", aspect = 1,
      ylab = "Difference in zinc concentrations (mg/L)"))


###################################################
### code chunk number 93: xmp0908
###################################################
str(xmp09.08)


###################################################
### code chunk number 94: xmp0908qqprt (eval = FALSE)
###################################################
## qqmath( ~ I(bottom-surface), xmp09.08)


###################################################
### code chunk number 95: xmp0908t
###################################################
with(xmp09.08, t.test(bottom,surface,paired=TRUE))


###################################################
### code chunk number 96: propdiffqq
###################################################
print(qqmath( ~ Difference, xmp09.09, type = c("g","p"),
      xlab = "Standard normal quantiles", aspect = 1,
      ylab = "Difference in proportion of time (%)"))


###################################################
### code chunk number 97: xmp0909
###################################################
str(xmp09.09)


###################################################
### code chunk number 98: xmp0909qqprt (eval = FALSE)
###################################################
## qqmath(~ Difference, xmp09.09)


###################################################
### code chunk number 99: xmp0909t
###################################################
with(xmp09.09, t.test(Before, After, paired = TRUE))


###################################################
### code chunk number 100: xmp0910
###################################################
Difference <- c(5,19,25,10,10,10,28,46,25,38,14,23,14)
qqnorm(Difference)
qqline(Difference)
t.test(Difference)


###################################################
### code chunk number 101: <xmp0910a
###################################################
with(xmp09.10, t.test(slide, digital, paired=TRUE))


###################################################
### code chunk number 102: <xmp1101a
###################################################
xmp11.01$brand = as.factor(xmp11.01$brand)
xmp11.01$treatment = as.factor(xmp11.01$treatment)


###################################################
### code chunk number 103: xmp1101
###################################################
with(xmp11.01, interaction.plot(treatment, brand, strength, col=2:4, lty=1))
with(xmp11.01, interaction.plot(treatment, brand, strength, col=2:4, lty=1))
with(xmp11.01, interaction.plot(brand, treatment, strength, col=2:4, lty=1))
with(xmp11.01, interaction.plot(brand, treatment, strength, col=2:4, lty=1))

anova(fm1 <- aov(strength ~ treatment + brand, xmp11.01))

TukeyHSD(fm1, which = "treatment") 

plot(fm1, which = 1:2)    


###################################################
### code chunk number 104: xmp1105
###################################################
str(xmp11.05)
with(xmp11.05, interaction.plot(humid, brand, power, col = 2:6, lty = 1))
anova(fm1 <- aov(power ~ humid + brand, data = xmp11.05))
TukeyHSD(fm1, which = "brand")
plot(TukeyHSD(fm1, which = "brand"))
plot(TukeyHSD(fm1, which = "brand"))


###################################################
### code chunk number 105: xmp1106
###################################################
str(xmp11.06)
anova(fm1 <- aov(Resp ~ Subject + Stimulus, data = xmp11.06))
with(xmp11.06, interaction.plot(Stimulus, Subject, Resp, col = 2:6, lty = 1))
TukeyHSD(fm1, which = "Stimulus")
plot(fm1,which = 1)

range(xmp11.06$Resp)
anova(fm2 <- aov(log(Resp) ~ Subject + Stimulus, data = xmp11.06))
with(xmp11.06, interaction.plot(Stimulus, Subject, log(Resp), col = 2:6, lty = 1))
plot(fm2,which = 1)
opar <- par(pty = 's')
plot(fm2,which = 2)
par(opar)
plot(TukeyHSD(fm2, which = "Stimulus"))
with(xmp11.06,interaction.plot(Stimulus, Subject, fitted(fm2), col=2:6, lty=1))


###################################################
### code chunk number 106: xmp1107
###################################################
str(xmp11.07)
xtabs(~ Variety + Density, data = xmp11.07)
xtabs(Yield ~ Variety + Density, data = xmp11.07)

xtabs(Yield ~ Variety + Density, xmp11.07)/ xtabs( ~ Variety + Density, xmp11.07)
anova(fm1 <- aov(Yield ~ Density * Variety, data = xmp11.07))
with(xmp11.07, interaction.plot(Density, Variety, Yield, col=2:6, lty=1))
anova(fm2 <- aov(Yield ~ Density + Variety, xmp11.07))
TukeyHSD(fm2)
model.tables(fm2, type = "means")


###################################################
### code chunk number 107: xmp1110
###################################################
str(xmp11.10)
xtabs(~ Period + Coat + Strain, xmp11.10)
anova(fm1 <- aov(Tempr ~ Period * Coat * Strain, xmp11.10))
anova(fm2 <- update(fm1, . ~ . - Period:Strain -  Period:Coat:Strain))


###################################################
### code chunk number 108: xmp1111
###################################################
str(xmp11.11)
xtabs(as.integer(humidity) ~ row + column, xmp11.11)
xtabs(abrasion ~ row + column, xmp11.11)
anova(fm1 <- aov(abrasion ~ row + column + humidity, xmp11.11))
model.tables(fm1, cterms = "humidity", type = "mean")
TukeyHSD(fm1, which = "humidity")


###################################################
### code chunk number 109: tpt
###################################################
2+2


###################################################
### code chunk number 110: showquit (eval = FALSE)
###################################################
## q()


###################################################
### code chunk number 111: datacommand (eval = FALSE)
###################################################
## data()


###################################################
### code chunk number 112: help (eval = FALSE)
###################################################
## help(pressure)
## ?pressure


###################################################
### code chunk number 113: installD6 (eval = FALSE)
###################################################
## install.packages('Devore7')


###################################################
### code chunk number 114: libraryD6
###################################################
library(Devore7)


###################################################
### code chunk number 115: DSlist (eval = FALSE)
###################################################
## data(package = 'Devore7')


###################################################
### code chunk number 116: struse
###################################################
data(xmp01.02)
str(xmp01.02)


###################################################
### code chunk number 117: strcat
###################################################
data(ex01.29)
str(ex01.29)


###################################################
### code chunk number 118: summarydemo
###################################################
summary(xmp01.02)
summary(ex01.29)


###################################################
### code chunk number 119: stemtemp
###################################################
data(xmp01.01)
str(xmp01.01)
stem(xmp01.01$temp)


###################################################
### code chunk number 120: stemtemp2
###################################################
attach(xmp01.01)
stem(temp)


###################################################
### code chunk number 121: vara
###################################################
a = c(1,2,3,4)
b = c(2,3,5,7)


###################################################
### code chunk number 122: dfab
###################################################
df = data.frame(a,b)
df


###################################################
### code chunk number 123: dfstacked
###################################################
stack(df)


