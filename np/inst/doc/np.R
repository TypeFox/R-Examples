### R code from vignette source 'np.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: np.Rnw:90-91
###################################################
options(prompt = "R> ", np.messages = FALSE, digits = 3)


###################################################
### code chunk number 2: np.Rnw:554-558
###################################################
library("np")
data("cps71")
model.par <- lm(logwage ~ age + I(age^2), data = cps71)
summary(model.par)


###################################################
### code chunk number 3: np.Rnw:570-576
###################################################
model.np <- npreg(logwage ~ age, 
                  regtype = "ll", 
                  bwmethod = "cv.aic",
                  gradients = TRUE,
                  data = cps71)
summary(model.np)


###################################################
### code chunk number 4: np.Rnw:592-593
###################################################
npsigtest(model.np)


###################################################
### code chunk number 5: np.Rnw:651-661 (eval = FALSE)
###################################################
## plot(cps71$age, cps71$logwage, xlab = "age", ylab = "log(wage)", cex=.1)
## lines(cps71$age, fitted(model.np), lty = 1, col = "blue")
## lines(cps71$age, fitted(model.par), lty = 2, col = " red")
## plot(model.np, plot.errors.method = "asymptotic")
## 
## plot(model.np, gradients = TRUE)
## lines(cps71$age, coef(model.par)[2]+2*cps71$age*coef(model.par)[3],
##       lty = 2,
##       col = "red")
## plot(model.np, gradients = TRUE, plot.errors.method = "asymptotic")


###################################################
### code chunk number 6: np.Rnw:666-691
###################################################
getOption("SweaveHooks")[["multifig"]]()
options(SweaveHooks = list(multifig = function() par(mfrow=c(2,2),mar=c(4,4,3,2))))

# Plot 1

plot(cps71$age,cps71$logwage,xlab="age",ylab="log(wage)",cex=.1)
lines(cps71$age,fitted(model.np),lty=1,col="blue")
lines(cps71$age,fitted(model.par),lty=2,col="red")

# Plot 2

plot(cps71$age,gradients(model.np),xlab="age",ylab="gradient",type="l",lty=1,col="blue")
lines(cps71$age,coef(model.par)[2]+2*cps71$age*coef(model.par)[3],lty=2,col="red")

# Plot 3

plot(cps71$age,fitted(model.np),xlab="age",ylab="log(wage)",ylim=c(min(fitted(model.np)-2*model.np$merr),max(fitted(model.np)+2*model.np$merr)),type="l")
lines(cps71$age,fitted(model.np)+2*model.np$merr,lty=2,col="red")
lines(cps71$age,fitted(model.np)-2*model.np$merr,lty=2,col="red")

# Plot 4

plot(cps71$age,gradients(model.np),xlab="age",ylab="gradient",ylim=c(min(gradients(model.np)-2*model.np$gerr),max(gradients(model.np)+2*model.np$gerr)),type="l",lty=1,col="blue")
lines(cps71$age,gradients(model.np)+2*model.np$gerr,lty=2,col="red")
lines(cps71$age,gradients(model.np)-2*model.np$gerr,lty=2,col="red")



###################################################
### code chunk number 7: np.Rnw:752-755
###################################################
cps.eval <- data.frame(age = seq(10,70, by=10))
predict(model.par, newdata = cps.eval)
predict(model.np, newdata = cps.eval)


###################################################
### code chunk number 8: np.Rnw:797-805
###################################################
data("wage1")
model.ols <- lm(lwage ~ female +
                married +
                educ +
                exper +
                tenure,
                data = wage1)
summary(model.ols)


###################################################
### code chunk number 9: np.Rnw:828-848
###################################################
model.ols <- lm(lwage ~ female +
                married +
                educ +
                exper +
                tenure,
                x = TRUE,
                y = TRUE,
                data = wage1)
X <- data.frame(wage1$female,
                wage1$married,
                wage1$educ,
                wage1$exper,
                wage1$tenure)
output <- npcmstest(model = model.ols, 
                    xdat = X, 
                    ydat = wage1$lwage, 
                    nmulti = 1,
                    tol = 0.1, 
                    ftol = 0.1)
summary(output)


###################################################
### code chunk number 10: np.Rnw:903-913
###################################################
#bw.all <- npregbw(formula = lwage ~ female +
#                  married +
#                  educ +
#                  exper +
#                  tenure,
#                  regtype = "ll",
#                  bwmethod = "cv.aic",
#                  data = wage1)
model.np <- npreg(bws = bw.all)
summary(model.np)


###################################################
### code chunk number 11: np.Rnw:946-981
###################################################
set.seed(123)
ii <- sample(seq(1, nrow(wage1)), replace=FALSE)
wage1.train <- wage1[ii[1:400],]
wage1.eval <- wage1[ii[401:nrow(wage1)],]
model.ols <- lm(lwage ~ female +
                married +
                educ +
                exper +
                tenure,
                data = wage1.train)
fit.ols <- predict(model.ols,
                   data = wage1.train,
                   newdata = wage1.eval)
pse.ols <- mean((wage1.eval$lwage - fit.ols)^2)
#bw.subset <- npregbw(formula = lwage ~ female +
#                     married +
#                     educ +
#                     exper +
#                     tenure,
#                     regtype = "ll",
#                     bwmethod = "cv.aic",
#                     data = wage1.train)
model.np <- npreg(bws = bw.subset)
fit.np <- predict(model.np,
                  data = wage1.train,
                  newdata = wage1.eval)
pse.np <- mean((wage1.eval$lwage - fit.np)^2)
bw.freq <- bw.subset
bw.freq$bw[1] <- 0
bw.freq$bw[2] <- 0
model.np.freq <- npreg(bws = bw.freq)
fit.np.freq <- predict(model.np.freq,
                       data = wage1.train,
                       newdata = wage1.eval)
pse.np.freq <- mean((wage1.eval$lwage - fit.np.freq)^2)


###################################################
### code chunk number 12: np.Rnw:1017-1020 (eval = FALSE)
###################################################
## plot(model.np,
##      plot.errors.method = "bootstrap",
##      plot.errors.boot.num = 25)


###################################################
### code chunk number 13: np.Rnw:1025-1028
###################################################
plot(model.np,
     plot.errors.method = "bootstrap",
     plot.errors.boot.num = 25)


###################################################
### code chunk number 14: np.Rnw:1073-1104
###################################################
data("birthwt", package = "MASS")
birthwt$low <- factor(birthwt$low)
birthwt$smoke <- factor(birthwt$smoke)
birthwt$race <- factor(birthwt$race)
birthwt$ht <- factor(birthwt$ht)
birthwt$ui <- factor(birthwt$ui)
birthwt$ftv <- factor(birthwt$ftv)
model.logit <- glm(low ~ smoke +
                   race +
                   ht +
                   ui +
                   ftv +
                   age +
                   lwt, 
                   family = binomial(link = logit),
                   data = birthwt)
model.np <- npconmode(low ~
                      smoke +
                      race +
                      ht +
                      ui +
                      ftv +
                      age +
                      lwt,
                      tol = 0.1,
                      ftol = 0.1,
                      data = birthwt)
cm <- table(birthwt$low, 
            ifelse(fitted(model.logit) > 0.5, 1, 0))
cm
summary(model.np)


###################################################
### code chunk number 15: np.Rnw:1132-1137
###################################################
data("faithful", package = "datasets")
f.faithful <- npudens(~ eruptions + waiting, data = faithful)
F.faithful <- npudist(~ eruptions + waiting, data = faithful)
summary(f.faithful)
summary(F.faithful)


###################################################
### code chunk number 16: np.Rnw:1142-1144 (eval = FALSE)
###################################################
## plot(f.faithful, xtrim = -0.2, view = "fixed", main = "")
## plot(F.faithful, xtrim = -0.2, view = "fixed", main = "")


###################################################
### code chunk number 17: np.Rnw:1149-1163
###################################################
getOption("SweaveHooks")[["multifig"]]()
options(SweaveHooks = list(multifig = function() par(mfrow=c(1,2),mar=rep(1.5,4))))
plot.output <- plot(f.faithful, xtrim=-0.2, view="fixed", main = "",plot.behavior="data")
# Retrieve data from plot() to create multiple figures
f <- matrix(plot.output$d1$dens, 50, 50)
plot.x1 <- unique(plot.output$d1$eval[,1])
plot.x2 <- unique(plot.output$d1$eval[,2])
persp(plot.x1,plot.x2,f,xlab="eruptions",ylab="waiting",zlab="Joint Density",col="lightblue", ticktype="detailed")

plot.output <- plot(F.faithful, xtrim = -0.2, view = "fixed", main="",plot.behavior="data")
# Retrieve data from plot() to create multiple figures
F <- matrix(plot.output$d1$dist, 50, 50)
plot.x1 <- unique(plot.output$d1$eval[,1])
plot.x2 <- unique(plot.output$d1$eval[,2])
persp(plot.x1,plot.x2,F,xlab="eruptions",ylab="waiting",zlab="Joint Distribution",col="lightblue", ticktype="detailed")


###################################################
### code chunk number 18: np.Rnw:1189-1200
###################################################
data("Italy")
fhat <- npcdens(gdp ~ year,
                tol = 0.1,
                ftol = 0.1,
                data = Italy)
summary(fhat)
Fhat <- npcdist(gdp ~ year,
                tol = 0.1,
                ftol = 0.1,
                data = Italy)
summary(Fhat)


###################################################
### code chunk number 19: np.Rnw:1207-1209 (eval = FALSE)
###################################################
## plot(fhat, view = "fixed", main = "", theta = 300, phi = 50)
## plot(Fhat, view = "fixed", main = "", theta = 300, phi = 50)


###################################################
### code chunk number 20: np.Rnw:1214-1229
###################################################
getOption("SweaveHooks")[["multifig"]]()
options(SweaveHooks = list(multifig = function() par(mfrow=c(1,2),mar=rep(1.25,4))))
plot.output <- plot(fhat, view="fixed", main="",plot.behavior="data")
# Retrieve data from plot() to create multiple figures
f <- matrix(plot.output$cd1$condens, 48, 50)
plot.y1 <- unique(plot.output$cd1$yeval)
plot.x1 <- unique(plot.output$cd1$xeval)
persp(as.integer(levels(plot.x1)),plot.y1,f,xlab="year",ylab="gdp",zlab="Conditional Density",col="lightblue", ticktype="detailed",theta=300,phi=50)

plot.output <- plot(Fhat, view="fixed", main="",plot.behavior="data")
# Retrieve data from plot() to create multiple figures
F <- matrix(plot.output$cd1$condist, 48, 50)
plot.y1 <- unique(plot.output$cd1$yeval)
plot.x1 <- unique(plot.output$cd1$xeval)
persp(as.numeric(plot.x1)+1951,plot.y1,F,xlab="year",ylab="gdp",zlab="Conditional Distribution",col="lightblue", ticktype="detailed",theta=300,phi=50)



###################################################
### code chunk number 21: np.Rnw:1259-1266
###################################################
bw <- npcdistbw(formula = gdp ~ year,
                tol = 0.1,
                ftol = 0.1,
                data = Italy)
model.q0.25 <- npqreg(bws = bw, tau = 0.25)
model.q0.50 <- npqreg(bws = bw, tau = 0.50)
model.q0.75 <- npqreg(bws = bw, tau = 0.75)


###################################################
### code chunk number 22: np.Rnw:1273-1280 (eval = FALSE)
###################################################
## plot(Italy$year, Italy$gdp, main = "", 
##      xlab = "Year", ylab = "GDP Quantiles")
## lines(Italy$year, model.q0.25$quantile, col = "red", lty = 1, lwd = 2)
## lines(Italy$year, model.q0.50$quantile, col = "blue", lty = 2, lwd = 2)
## lines(Italy$year, model.q0.75$quantile, col = "red", lty = 3, lwd = 2)
## legend(ordered(1951), 32, c("tau = 0.25", "tau = 0.50", "tau = 0.75"), 
##        lty = c(1, 2, 3), col = c("red", "blue", "red"))


###################################################
### code chunk number 23: np.Rnw:1289-1298
###################################################
plot(Italy$year, Italy$gdp, 
     main = "",
     xlab = "Year", 
     ylab = "GDP Quantiles")
lines(Italy$year, model.q0.25$quantile, col = "red", lty = 1, lwd = 2)
lines(Italy$year, model.q0.50$quantile, col = "blue", lty = 2, lwd = 2)
lines(Italy$year, model.q0.75$quantile, col = "red", lty = 3, lwd = 2)
legend(ordered(1951), 32, c("tau = 0.25", "tau = 0.50", "tau = 0.75"), 
       lty = c(1, 2, 3), col = c("red", "blue", "red"))


###################################################
### code chunk number 24: np.Rnw:1341-1347
###################################################
model.pl <- npplreg(lwage ~ female +
                    married +
                    educ +
                    tenure | exper,
                    data = wage1)
summary(model.pl)


###################################################
### code chunk number 25: np.Rnw:1367-1379
###################################################
model.index <- npindex(low ~
                       smoke +
                       race +
                       ht +
                       ui +
                       ftv +
                       age +
                       lwt,
                       method = "kleinspady",
                       gradients = TRUE,
                       data = birthwt)
summary(model.index)


###################################################
### code chunk number 26: np.Rnw:1399-1407
###################################################
model <- npindex(lwage ~ female +
                 married +
                 educ +
                 exper +
                 tenure,
                 data = wage1,
                 nmulti = 1)
summary(model)


###################################################
### code chunk number 27: np.Rnw:1433-1452
###################################################
model.ols <- lm(lwage ~ female +
                married +
                educ +
                exper +
                tenure,
                data = wage1)
wage1.augmented <- wage1
wage1.augmented$dfemale <- as.integer(wage1$female == "Male")
wage1.augmented$dmarried <- as.integer(wage1$married == "Notmarried")
model.scoef <- npscoef(lwage ~ dfemale +
                       dmarried +
                       educ +
                       exper +
                       tenure | female,
                       betas = TRUE,
                       data = wage1.augmented)
summary(model.scoef)
colMeans(coef(model.scoef))
coef(model.ols)


###################################################
### code chunk number 28: np.Rnw:1478-1480
###################################################
fit.lc <- npksum(txdat = cps71$age, tydat = cps71$logwage, bws = 2)$ksum/
  npksum(txdat = cps71$age, bws = 2)$ksum


###################################################
### code chunk number 29: np.Rnw:1485-1487 (eval = FALSE)
###################################################
## plot(cps71$age, cps71$logwage, xlab = "Age", ylab = "log(wage)")
## lines(cps71$age, fit.lc, col = "blue")


###################################################
### code chunk number 30: np.Rnw:1492-1494
###################################################
plot(cps71$age, cps71$logwage, xlab = "Age", ylab = "log(wage)", cex=.1)
lines(cps71$age, fit.lc, col = "blue")


