### R code from vignette source 'vignette.Rnw'

###################################################
### code chunk number 1: vignette.Rnw:36-37 (eval = FALSE)
###################################################
## install.packages("Rearrangement")


###################################################
### code chunk number 2: vignette.Rnw:41-42
###################################################
library(Rearrangement)


###################################################
### code chunk number 3: vignette.Rnw:46-47 (eval = FALSE)
###################################################
## help(package="Rearrangement")


###################################################
### code chunk number 4: vignette.Rnw:51-52 (eval = FALSE)
###################################################
## help(rearrangement)


###################################################
### code chunk number 5: vignette.Rnw:74-76
###################################################
data(GrowthChart)
attach(GrowthChart)


###################################################
### code chunk number 6: vignette.Rnw:79-80
###################################################
library(splines)


###################################################
### code chunk number 7: vignette.Rnw:83-89
###################################################
ages <- unique(sort(age))
aknots <- c(3, 5, 8, 10, 11.5, 13, 14.5, 16, 18)
splines_age <- bs(age,kn=aknots)
sformula <- height~splines_age
sfunc <- approxfun(age,lm(sformula)$fitted.values)
splreg <- sfunc(ages)


###################################################
### code chunk number 8: vignette.Rnw:92-93
###################################################
rsplreg <- rearrangement(list(ages),splreg)


###################################################
### code chunk number 9: vignette.Rnw:96-100
###################################################
plot(age,height,pch=21,bg='gray',cex=.5,xlab="Age (years)",ylab="Height (cms)",main="CEF (Regression Splines)",col='gray')
lines(ages,splreg,col='red',lwd=3)
lines(ages,rsplreg,col='blue',lwd=2)
legend("topleft",c('Original','Rearranged'),lty=1,col=c('red','blue'),bty='n')


###################################################
### code chunk number 10: vignette.Rnw:115-116 (eval = FALSE)
###################################################
## rearrangement(x, y, n = 1000, stochastic = FALSE, avg = TRUE, order = 1:length(x))


###################################################
### code chunk number 11: vignette.Rnw:130-131
###################################################
f <- function(x1,x2){(x1 - 1)*(x1 - 2)*(x1 - 3)*(x2 - 1)*(x2 - 2)*(x2 - 3)}


###################################################
### code chunk number 12: vignette.Rnw:134-140
###################################################
x1 <- seq(0, 5, by = .5)
x2 <- x1
l <- NULL
for (i in x2){l <- c(l, f(x1, i))}
y <- matrix(l, nrow = length(x1))
x <- list(x1,x2)


###################################################
### code chunk number 13: vignette.Rnw:143-144
###################################################
ry <- rearrangement(x,y)


###################################################
### code chunk number 14: vignette.Rnw:148-154
###################################################
par(mfrow=c(1,2))
persp(x1,x2,y,col='lightgreen',theta=315,phi=25,r=6,shade=.5,font.lab=15,cex.lab=3)
title("Before Rearrangement",cex.main=3,font.main=30)

persp(x1,x2,ry,col='lightgreen',theta=315,phi=25,r=6,shade=.5,font.lab=15,cex.lab=3)
title("After Rearrangement",cex.main=3,font.main=30)


###################################################
### code chunk number 15: vignette.Rnw:156-157
###################################################
par(mfrow=c(1,1))


###################################################
### code chunk number 16: vignette.Rnw:197-199
###################################################
data(GrowthChart)
attach(GrowthChart)


###################################################
### code chunk number 17: vignette.Rnw:203-209
###################################################
##Normalize the ages to the interval [0, 2*pi]
nage <- 2 * pi * (age - min(age)) / (max(age) - min(age))
formula <- height ~ I(sin(nage)) + I(cos(nage)) + I(sin(2*nage)) + I(cos(2*nage)) + I(sin(3*nage)) + I(cos(3*nage))  + I(sin(4*nage)) + I(cos(4*nage))
j <- simconboot(nage, height, lm, formula)
class(j)
names(j)


###################################################
### code chunk number 18: vignette.Rnw:213-215
###################################################
plot(j, border=NA, col='darkgray', xlab='Age (years)', ylab='Height (cms)', xaxt = "n")
axis(1, at = seq(-2*pi*min(age)/(max(age)-min(age)), 2*pi+1, by=5*2*pi/(max(age)-min(age))), label = seq(0, max(age)+1, by=5))


###################################################
### code chunk number 19: vignette.Rnw:228-229
###################################################
k <- rconint(j)


###################################################
### code chunk number 20: vignette.Rnw:232-248
###################################################
ages <- unique(sort(age))
ffunc <- approxfun(age,lm(formula)$fitted.values)
freg <- ffunc(ages)
rfreg <- rearrangement(list(ages),freg)

plot(k, border=NA, col='#2A2A2A', xlab='Age (years)', ylab='Height (cms)', xaxt = "n")
axis(1, at = seq(-2*pi*min(age)/(max(age)-min(age)), 2*pi+1, by=5*2*pi/(max(age)-min(age))), label = seq(0, max(age)+1, by=5))
polygon.conint(j, border=NA, col='#D2D2D2')
polygon.conint(k, border=NA, col='#2A2A2A', density=50)
points(nage,height, cex = .5, col = '#7E7E7E')

nages <- unique(sort(nage))
lines(nages,freg, col='red',lwd=2)
lines(nages,rfreg, col='blue',lwd=2)
legend("topleft",c('95% CI Original','95% CI Rearranged'),lty=1,lwd=7,col=c('#D2D2D2','#2A2A2A'),bty='n')
legend("topleft",c('95% CI Original','95% CI Rearranged'),lty=1,col=c('red','blue'),bty='n')


###################################################
### code chunk number 21: vignette.Rnw:256-263
###################################################
library(quantreg)
ages <- unique(sort(age))

## Univariate
j2 <- simconboot(age, height, lprq2, formula=0,B=20,h=1,xx=ages,tau=0.5)
k2 <- rconint(j2)
rcqf50 <- rearrangement(data.frame(j2$sortedx),j2$cef)


###################################################
### code chunk number 22: vignette.Rnw:269-278
###################################################
plot(age,height,xlab="Age (years)",ylab = "Height(cms)",col='grey80')
polygon.conint(j2, border=NA, col='#D2D2D2')
polygon.conint(k2, border=NA, col='#2A2A2A', density=50)
lines(j2$sortedx,j2$cef,lty=1,lwd=2.5,col="tomato2")
lines(j2$sortedx,rcqf50,lty=1,lwd=1,col="blue")

legend("topleft",c('95% CI Original','95% CI Rearranged'),lty=1,lwd=7,col=c('#D2D2D2','#2A2A2A'),bty='n')
legend("topleft",c('95% CI Original','95% CI Rearranged'),lty=1,col=c('red','blue'),bty='n')
title(main="CQF (Local linear, h=1) ")


###################################################
### code chunk number 23: vignette.Rnw:298-406 (eval = FALSE)
###################################################
## 
## formula0	<- height ~ age + I((age >= 5)*(age - 5)) + I((age >= 10)*(age - 10))  + I((age >= 15)*(age - 15));
## cef		<- approxfun(age, lm(formula0)$fitted.values);
## ages		<- unique(sort(age));
## mheight 	<- cef(ages);
## 
## ###### Compute Locally Linear Mean Regression ###########################################
## # (Modification of Koenker's lprq with box kernel)
## 
## # Residuals
## residuals0	<- lm(formula0)$residuals;
## 
## set.seed(8);   
## n		<- length(height); 	# sample size
## B		<- 20; 		# Number of replications
## h		<- .5
## 
## formula1 	<- heightb ~ as.factor(age);
## 
## nknots 		<- 9;
## knots_age	<- c(3.0, 5.0, 8.0, 10, 11.5, 13, 14.5, 16, 18);
## splines_age     <- bs(age, degree = 3, intercept = FALSE, knots = knots_age);
## formula2	<- heightb ~ splines_age;
## 
## nage		<- 2 * pi * (age - min(age)) / (max(age) - min(age));
## formula5	<- heightb ~ nage + I(nage^2) + I(sin(nage)) + I(cos(nage)) + I(sin(2*nage)) + I(cos(2*nage));
## 
## mheight0	<- NULL;
## 
## mheight2	<- NULL;
## rmheight2	<- NULL;
## 
## mheight5	<- NULL;
## rmheight5	<- NULL;
## 
## mheight7	<- NULL;
## rmheight7	<- NULL;
## 
## mheight8	<- NULL;
## rmheight8	<- NULL;
## 
## for (s in 1:B) {;
## 
## res 	<- rnorm(n, mean = mean(residuals0), sd = sd(residuals0));
## heightb <- cef(age) + res;
## 
## mheight0	<- rbind(mheight0, mheight);
## 
## 
## cef2		<- approxfun(age, lm(formula2)$fitted.values);
## 	mheight2	<- rbind(mheight2, cef2(ages));
## rmheight2	<- rbind(rmheight2, rearrangement(list(ages), mheight2[s, ]));
## 
## cef5		<- approxfun(age, lm(formula5)$fitted.values);
## 	mheight5	<- rbind(mheight5, cef5(ages));
## rmheight5	<- rbind(rmheight5, rearrangement(list(ages), mheight5[s, ]));
## 
## cef7		<- lplm(age, heightb, h, ages)$fitted.values;
## 	mheight7	<- rbind(mheight7, cef7);
## rmheight7	<- rbind(rmheight7, rearrangement(list(ages), mheight7[s, ]));
## 
## 	cef8		<- lclm(age, heightb, h, ages)$fitted.values;
## 	mheight8	<- rbind(mheight8, cef8);
## 	rmheight8	<- rbind(rmheight8, rearrangement( list(ages), mheight8[s, ]));
## 
## }
## 
## table <- matrix(0, nrow = 3, ncol = 12, dimnames = list(c("L1", "L2", "Linf"), c("Splines", "Rearranged", "Ratio (R/O)","Fourier", "Rearranged", "Ratio (R/O)","Local Poly.", "Rearranged", "Ratio (R/O)","Kernel (h=1)", "Rearranged", "Ratio (R/O)")));
## 
## 
## table[1,1] <- mean(apply(abs(mheight2 - mheight0), 1, mean));
## table[1,2] <- mean(apply(abs(rmheight2 - mheight0), 1, mean));
## table[1,3] <- table[1,2]/table[1,1];
## table[1,4] <- mean(apply(abs(mheight5 - mheight0), 1, mean));
## table[1,5] <- mean(apply(abs(rmheight5 - mheight0), 1, mean));
## table[1,6] <- table[1,5]/table[1,4];
## table[1,7] <- mean(apply(abs(mheight7 - mheight0), 1, mean));
## table[1,8] <- mean(apply(abs(rmheight7 - mheight0), 1, mean));
## table[1,9] <- table[1,8]/table[1,7];
## table[1,10] <- mean(apply(abs(mheight8 - mheight0), 1, mean));
## table[1,11] <- mean(apply(abs(rmheight8 - mheight0), 1, mean));
## table[1,12] <- table[1,11]/table[1,10];
## 
## table[2,1] <- mean(apply(abs(mheight2 - mheight0)^2, 1, mean)^(1/2));
## table[2,2] <- mean(apply(abs(rmheight2 - mheight0)^2, 1, mean)^(1/2));
## table[2,3] <- table[2,2]/table[2,1];
## table[2,4] <- mean(apply(abs(mheight5 - mheight0)^2, 1, mean)^(1/2));
## table[2,5] <- mean(apply(abs(rmheight5 - mheight0)^2, 1, mean)^(1/2));
## table[2,6] <- table[2,5]/table[2,4];
## table[2,7] <- mean(apply(abs(mheight7 - mheight0)^2, 1, mean)^(1/2));
## table[2,8] <- mean(apply(abs(rmheight7 - mheight0)^2, 1, mean)^(1/2));
## table[2,9] <- table[2,8]/table[2,7];
## table[2,10] <- mean(apply(abs(mheight8 - mheight0)^2, 1, mean)^(1/2));
## table[2,11] <- mean(apply(abs(rmheight8 - mheight0)^2, 1, mean)^(1/2));
## table[2,12] <- table[2,11]/table[2,10];
## 
## table[3,1] <- mean(apply(abs(mheight2 - mheight0), 1, max));
## table[3,2] <- mean(apply(abs(rmheight2 - mheight0), 1, max));
## table[3,3] <- table[3,2]/table[3,1];
## table[3,4] <- mean(apply(abs(mheight5 - mheight0), 1, max));
## table[3,5] <- mean(apply(abs(rmheight5 - mheight0), 1, max));
## table[3,6] <- table[3,5]/table[3,4];
## table[3,7] <- mean(apply(abs(mheight7 - mheight0), 1, max));
## table[3,8] <- mean(apply(abs(rmheight7 - mheight0), 1, max));
## table[3,9] <- table[3,8]/table[3,7];
## table[3,10] <- mean(apply(abs(mheight8 - mheight0), 1, max));
## table[3,11] <- mean(apply(abs(rmheight8 - mheight0), 1, max));
## table[3,12] <- table[3,11]/table[3,10];


###################################################
### code chunk number 24: vignette.Rnw:409-410 (eval = FALSE)
###################################################
## table


