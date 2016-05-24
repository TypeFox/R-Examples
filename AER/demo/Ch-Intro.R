###################################################
### chunk number 1: setup
###################################################
options(prompt = "R> ", continue = "+  ", width = 64,
  digits = 4, show.signif.stars = FALSE, useFancyQuotes = FALSE)

options(SweaveHooks = list(onefig =   function() {par(mfrow = c(1,1))},
                           twofig =   function() {par(mfrow = c(1,2))},                           
                           threefig = function() {par(mfrow = c(1,3))},
                           fourfig =  function() {par(mfrow = c(2,2))},
			   sixfig =   function() {par(mfrow = c(3,2))}))

library("AER")

set.seed(1071)


###################################################
### chunk number 2: journals-data
###################################################
data("Journals", package = "AER")


###################################################
### chunk number 3: journals-dim
###################################################
dim(Journals)
names(Journals)


###################################################
### chunk number 4: journals-plot eval=FALSE
###################################################
## plot(log(subs) ~ log(price/citations), data = Journals)


###################################################
### chunk number 5: journals-lm eval=FALSE
###################################################
## j_lm <- lm(log(subs) ~ log(price/citations), data = Journals)
## abline(j_lm)


###################################################
### chunk number 6: journals-lmplot
###################################################
plot(log(subs) ~ log(price/citations), data = Journals)
j_lm <- lm(log(subs) ~ log(price/citations), data = Journals)
abline(j_lm)


###################################################
### chunk number 7: journals-lm-summary
###################################################
summary(j_lm)


###################################################
### chunk number 8: cps-data
###################################################
data("CPS1985", package = "AER")
cps <- CPS1985


###################################################
### chunk number 9: cps-data1 eval=FALSE
###################################################
## data("CPS1985", package = "AER")
## cps <- CPS1985


###################################################
### chunk number 10: cps-reg
###################################################
library("quantreg")
cps_lm <- lm(log(wage) ~ experience + I(experience^2) +
  education, data = cps)
cps_rq <- rq(log(wage) ~ experience + I(experience^2) +
  education, data = cps, tau = seq(0.2, 0.8, by = 0.15))


###################################################
### chunk number 11: cps-predict
###################################################
cps2 <- data.frame(education = mean(cps$education),
  experience = min(cps$experience):max(cps$experience))
cps2 <- cbind(cps2, predict(cps_lm, newdata = cps2,
  interval = "prediction"))
cps2 <- cbind(cps2,
  predict(cps_rq, newdata = cps2, type = ""))  


###################################################
### chunk number 12: rq-plot eval=FALSE
###################################################
## plot(log(wage) ~ experience, data = cps)
## for(i in 6:10) lines(cps2[,i] ~ experience,
##   data = cps2, col = "red")


###################################################
### chunk number 13: rq-plot1
###################################################
plot(log(wage) ~ experience, data = cps)
for(i in 6:10) lines(cps2[,i] ~ experience,
  data = cps2, col = "red")


###################################################
### chunk number 14: srq-plot eval=FALSE
###################################################
## plot(summary(cps_rq))


###################################################
### chunk number 15: srq-plot1
###################################################
plot(summary(cps_rq))


###################################################
### chunk number 16: bkde-fit
###################################################
library("KernSmooth")
cps_bkde <- bkde2D(cbind(cps$experience, log(cps$wage)),
  bandwidth = c(3.5, 0.5), gridsize = c(200, 200))


###################################################
### chunk number 17: bkde-plot eval=FALSE
###################################################
## image(cps_bkde$x1, cps_bkde$x2, cps_bkde$fhat, 
##   col = rev(gray.colors(10, gamma = 1)),
##   xlab = "experience", ylab = "log(wage)")
## box()
## lines(fit ~ experience, data = cps2)
## lines(lwr ~ experience, data = cps2, lty = 2)
## lines(upr ~ experience, data = cps2, lty = 2)


###################################################
### chunk number 18: bkde-plot1
###################################################
image(cps_bkde$x1, cps_bkde$x2, cps_bkde$fhat, 
  col = rev(gray.colors(10, gamma = 1)),
  xlab = "experience", ylab = "log(wage)")
box()
lines(fit ~ experience, data = cps2)
lines(lwr ~ experience, data = cps2, lty = 2)
lines(upr ~ experience, data = cps2, lty = 2)


###################################################
### chunk number 19: install eval=FALSE
###################################################
## install.packages("AER")


###################################################
### chunk number 20: library
###################################################
library("AER")


###################################################
### chunk number 21: objects
###################################################
objects()


###################################################
### chunk number 22: search
###################################################
search()


###################################################
### chunk number 23: assignment
###################################################
x <- 2
objects()


###################################################
### chunk number 24: remove
###################################################
remove(x)
objects()


###################################################
### chunk number 25: log eval=FALSE
###################################################
## log(16, 2)
## log(x = 16, 2)
## log(16, base = 2)
## log(base = 2, x = 16)


###################################################
### chunk number 26: q eval=FALSE
###################################################
## q()


###################################################
### chunk number 27: apropos
###################################################
apropos("help")


