
## @knitr echo=FALSE, out.width=singlewide
## make a simple simulation of y = 1 + 2 x + e
## mark regression line, and actual relationship
x <- rep(1:10,5)
y <- rnorm(length(x), 1 + 2*x, 3)

plot(y ~ x, xlim=c(0,10), bty="l")
res <- lm(y~x)
abline(a=1,b=2,lty=1)
abline(res,lty=2)
points(mean(x),mean(y),pch=15,col=gray(.8),cex=3)
legend(0,.95*max(y),
       legend=c("Actual model","Estimated regression line"),
       lty=1:2,
       y.intersp = 2, bty="n")


## @knitr 
res.mhr <- lm(maxrate ~ age, data=heartrate)
res.mhr


## @knitr 
208.36 - 0.76 * 50


## @knitr 
n <- length(heartrate$age)
sum( resid(res)^2 ) / (n-2) 


## @knitr 
deviance(res)/ (n - 2) 


## @knitr 
library(MASS)                 # loads data set


## @knitr 
plot(MPG.highway ~ Horsepower, data = Cars93)
res <- lm(MPG.highway ~ Horsepower, data = Cars93)
res
predict(res, newdata=data.frame(Horsepower=225))


## @knitr 
f <- MPG.highway ~ Weight
plot(f, data=Cars93)
res <- lm(f, data=Cars93)
res
predict(res, newdata=data.frame(Weight=c(2524, 6400)))


## @knitr 
f <- Max.Price ~ Min.Price
plot(f, data=Cars93)
res <- lm(f, data=Cars93)
abline(res)
res


## @knitr 
rlm(f, data=Cars93)


## @knitr 
cars <- Cars93[,sapply(Cars93, function(x) !is.factor(x))]
pairs(cars)


## @knitr 
lm(attendance ~ wins,MLBattend)
27345*10


## @knitr 
age.two <- c(39, 30, 32, 34, 35, 36, 36, 30)
adult <- c(71, 63, 63, 67, 68, 68, 70, 64)
lm(adult ~ age.two)
35.179 + .929*33


## @knitr 
plot(child ~ parent, data=galton)
res <- lm(child ~ parent)
abline(res)
points(jitter(parent,5), jitter(child,5), pch=3)
plot(parent, child)
res <- lm(child ~ parent)
abline(res)
points(jitter(parent,5), jitter(child,5), pch=3)
res
cor(parent,child)


## @knitr 
lm(scale(child) ~ scale(parent), data=galton)
with(galton, cor(scale(child), scale(parent)))


## @knitr echo=FALSE, out.width=doublewide
x <- rep(1:10, 10)
y <- rnorm(length(x), x, 1)
plot(x, y, pch=16)

betahat0 <- numeric(); betahat1 <- numeric()
for(i in 1:100) {
  y <- rnorm(length(x), x, 1)
  res <- lm(y ~ x)
  betahat0[i] = coef(res)[1]
  betahat1[i] = coef(res)[2]
  abline(res, col=rgb(.75, .75, .75,.5)) #color?
}
plot(betahat0, betahat1)
hist(betahat0)
hist(betahat1)


## @knitr 
res.mhr <- lm(maxrate ~ age, data=heartrate)

betahat0 <- coef(res.mhr)[1]       # first coefficient
n <- nrow(heartrate)
sigmahat <- sqrt( sum(resid(res.mhr)^2) / (n - 2))
SE <- with(heartrate, 
           sigmahat*sqrt(sum(age^2) / (n*sum((age - mean(age))^2))) 
           )
SE
tstar <- qt(1 - 0.05/2, df=n - 2)

betahat0 + c(-1, 1) * tstar * SE


## @knitr cache=FALSE
summary(res.mhr)


## @knitr 
betahat1 <- -0.7595           # read from summary
SE <- 0.0561                  # read from summary
tstar <- qt(1 - 0.05/2, df=n - 2)
betahat1 + c(-1, 1) * tstar * SE


## @knitr cache=FALSE
coef(res.mhr)


## @knitr 
coef(summary(res.mhr))


## @knitr 
coef(summary(res.mhr))["age", "Std. Error"]


## @knitr 
mu0 <- -1
T.obs <- (betahat1 - mu0)/SE
T.obs
2*pt(abs(T.obs), df=n-2, lower.tail=FALSE)   


## @knitr 
sigma2 <- sum(resid(res.mhr)^2) / res.mhr$df.residual
sqrt(sigma2)                            # sigma hat


## @knitr 
(-0.7595 / 0.0561)^2


## @knitr 
anova(res.mhr)


## @knitr 
short_summary <- function(x) {
  x <- summary(x)
  cmat <- coef(x)
  printCoefmat(cmat)
}


## @knitr 
res.mhr <- lm(maxrate ~ age, data=heartrate)
predict(res.mhr, newdata=data.frame(age=42))


## @knitr echo=FALSE, out.width=doublewide
## upper left
x <- rep(1:10,4)
y <- rnorm(40, mean=5*sin(x), sd=1)
plot(y ~ x, main="Scatterplot")
abline(lm(y~x))
## upper right
x <- rep(1:10,4)
y <- rnorm(40,mean = x + .05*sin(x),sd=.01) # small trend
res <- lm(y~x)
plot(fitted(res),resid(res), main="fitted-resid")
## lower left
x <- rep(1:10,4)
y <- rnorm(40, mean = 1 + 1/2*x, sd = x/10)
res <- lm(y ~ x)
plot(fitted(res),resid(res), main="fitted-resid") 
## lower right
x <- rep(1:10, 4)
epsilon <- rnorm(40, mean=0, sd=1)
y <- 1 + 2*x + cumsum(epsilon) # cumsum() correlates errors
res <- lm(y ~ x)
tmp <- resid(res)
n <- length(tmp)
plot(tmp[-n], tmp[-1], main="lag plot")         # lag plot


## @knitr eval=FALSE
## x <- rep(1:10,4)
## y <- rnorm(40, mean=5*sin(x), sd=1)
## plot(y ~ x)
## abline(lm(y ~ x))


## @knitr eval=FALSE
## x <- rep(1:10, 4)
## y <- rnorm(40, mean=x + 0.05 * sin(x), sd=0.01) # small trend
## res <- lm(y ~ x)
## plot(fitted(res), resid(res))


## @knitr fig.keep="none"
x <- rep(1:10, 4)
y <- rnorm(40, mean=1 + 1/2*x, sd=x/10)
res <- lm(y ~ x)
plot(fitted(res), resid(res))


## @knitr eval=FALSE
## x <- rep(1:10, 4)
## epsilon <- rnorm(40, mean=0, sd=1)
## y <- 1 + 2*x + cumsum(epsilon) # cumsum() correlates errors
## res <- lm(y ~ x)
## tmp <- resid(res)
## n <- length(tmp)
## plot(tmp[-n], tmp[-1])         # lag plot


## @knitr emission-with-cooks, eval=FALSE
## res <- lm(CO2 ~ perCapita, emissions)
## plot(CO2 ~ perCapita, emissions,
##      cex=10*sqrt(cooks.distance(res)),
##      main=expression(                       # make subscript on C02
##        paste("bubble plot of ",CO[2],
##              " emissions by per capita GDP")
##        ))


## @knitr echo=FALSE, out.width=singlewide
res <- lm(CO2 ~ perCapita, emissions)
plot(CO2 ~ perCapita, emissions,
     cex=10*sqrt(cooks.distance(res)),
     main=expression(                       # make subscript on C02
       paste("bubble plot of ",CO[2],
             " emissions by per capita GDP")
       ))


## @knitr echo=FALSE, out.width=doublewide
res <- lm(maxrate ~ age, data=heartrate)
plot(res,which=1)
plot(res,which=2)
plot(res,which=3)
plot(res,which=5)


## @knitr 
pred.res <- predict(res.mhr, int = "pred")
head(pred.res, n=3)


## @knitr 
head(pred.res[, 2])                 # the 'lwr' column


## @knitr 
age.sort <- sort(unique(heartrate$age))
pred.res <- predict(res.mhr, newdata = data.frame(age = age.sort),
                    int="pred")
pred.res[,2]


## @knitr age-v-mhr-with-prediction, eval=FALSE
## plot(maxrate ~ age, data=heartrate)
## abline(res.mhr)
## lines(age.sort, pred.res[,2], lty=2) # lower curve
## lines(age.sort, pred.res[,3], lty=2) # upper curve


## @knitr echo=FALSE, out.width=singlewide
plot(maxrate ~ age, data=heartrate)
abline(res.mhr)
lines(age.sort, pred.res[,2], lty=2) # lower curve
lines(age.sort, pred.res[,3], lty=2) # upper curve


## @knitr 
price <- c(300, 250, 400, 550, 317, 389, 425, 289, 389, 559)
no.bed <- c(3, 3, 4, 5, 4, 3, 6, 3, 4, 5)
res <- lm(price ~ no.bed)
summary(res)
T <- (73.1 - 60)/23.8
pt(T, df=8, lower.tail=FALSE)


## @knitr echo=FALSE
tmp <- cbind(c(600, 1000, 1250, 1600, 1800, 2100, 2500, 2900),
             c(56, 54, 56, 50, 47, 49, 47, 45))
write.table(tmp, "tmp.txt")


## @knitr 
x <- read.table(file="tmp.txt")
names(x) <- c("elev", "Temp")
res <- lm(Temp ~ elev, x)
summary(res)
T <- (-0.005115 - (-0.00534)) / 0.000921
T
2*pt(T, df=6, lower.tail=FALSE) # two-sided


## @knitr 
no.beers <- c(5, 2, 9, 8, 3, 7, 3, 5, 3, 5)
BAL <- c(0.10, 0.03, 0.19, 0.12, 0.04, 0.095, 0.07, 0.06, 0.02, 0.05)
res <- lm(BAL ~ no.beers)
summary(res)
T <- (0.01920 - .02)/0.00351
2*pt(T, df = length(no.beers)-2)


## @knitr 
no.beers <- c(5, 2, 9, 8, 3, 7, 3, 5, 3, 5)
BAL <- c(0.10, 0.03, 0.19, 0.12, 0.04, 0.095, 0.07, 0.06, 0.02, 0.05)
short_summary(lm(BAL ~ no.beers))


## @knitr 
res <- lm(y2000 ~ y1970, data=homedata)
predict(res, newdata=dataframe(y1970=80000))
predict(res, newdata=data.frame(y1970=80000))
plot(resid(res))                  # simple plot
plot(homedata$y1970, resid(res))  # shows spread


## @knitr 
years <- 1952:1962
seal.counts <- c(724, 756, 920, 1392, 1392, 1448, 1212, 1672, 
  2068, 1980, 2116)
plot(seal.counts ~ years)
res <- lm(seal.counts ~ years)
predict(res, newdata=data.frame(years=1963))


## @knitr 
by.dist <- split(best.times, as.factor(best.times$Dist))    


## @knitr fig.keep="none"
plot(Time ~ age, by.dist[["800"]])    


## @knitr 
lm(scale(Time) ~ age, by.dist[["800"]], subset = age < 70)


## @knitr 
sapply(c("100", "400", "10000"), function(distance) {
  res <- lm(scale(Time) ~ age, by.dist[[distance]], subset = age < 70)
  coef(res)[2]
})


## @knitr 
res <- lm(child ~ parent, data=galton)
summary(res)
T <- (0.6463-1)/0.0411
2*pt(T, df=926)


## @knitr 
age.sort <- sort(unique(heartrate$age))
pred.res <- predict(res.mhr, newdata=data.frame(age=age.sort), int="conf")
plot(maxrate ~ age, heartrate); abline(res)
lines(age.sort,pred.res[,3], lty=2)
lines(age.sort,pred.res[,2], lty=2)


## @knitr 
res <- lm(log(lab.defect) ~ log(field.defect), data=alaska.pipeline)
plot(resid(res) ~ log(field.defect), data=alaska.pipeline)


## @knitr fig.keep="none"
m <- 200
x <- rep(1:10, 4)
res <- replicate(m, {
  y <- rnorm(40, 1 + 2*x, 3)
  coef(lm(y ~ x))
})
plot(res[1,], res[2,])


## @knitr 
library(ellipse)
res <- lm(Deflection ~ Load, data=deflection)
plot(ellipse(res), type="l")


## @knitr eval=FALSE
## t.test(temperature ~ factor(gender), data=normtemp)


## @knitr eval=FALSE
## lm(temperature ~ factor(gender), data=normtemp)


## @knitr 
x <- 1:10
y <- rchisq(10,3)
z <- 1 + x + y + rnorm(10)
lm(z ~ x + y)


## @knitr 
res.lm <- lm(wt ~ gestation + age + ht + wt1 + dage + dht + dwt , 
 data = babies,
 subset= gestation < 350 & age < 99 & ht < 99 & wt1 < 999 &
 dage < 99 & dht < 99 & dwt < 999)


## @knitr fig.keep="none"
plot(fitted(res.lm), resid(res.lm)) 


## @knitr 
init.h <- c(600,700,800,950,1100,1300,1500)
h.d <- c(253, 337, 395, 451, 495, 534, 573)
res.lm  <- lm(h.d ~ init.h)
res.lm2 <- update(res.lm,  . ~ . + I(init.h^2))
res.lm3 <- update(res.lm2, . ~ . + I(init.h^3))


## @knitr 
polynomial <- Vectorize(function(x, ps) {
  n <- length(ps)
  sum(ps*x^(1:n-1))
}, "x")


## @knitr galileo-fits, eval=FALSE
## plot(h.d ~ init.h)
## curve(polynomial(x, coef(res.lm )), add=TRUE, lty=1)
## curve(polynomial(x, coef(res.lm2)), add=TRUE, lty=2)
## curve(polynomial(x, coef(res.lm3)), add=TRUE, lty=3)
## legend(1200, 400, legend=c("linear", "quadratic", "cubic"), lty=1:3)


## @knitr echo=FALSE, out.width=singlewide, cache=FALSE
plot(h.d ~ init.h)
curve(polynomial(x, coef(res.lm )), add=TRUE, lty=1)
curve(polynomial(x, coef(res.lm2)), add=TRUE, lty=2)
curve(polynomial(x, coef(res.lm3)), add=TRUE, lty=3)
legend(1200, 400, legend=c("linear", "quadratic", "cubic"), lty=1:3)


## @knitr echo=FALSE
## Galileo data ...
init.h <- c(600,700,800,950,1100,1300,1500)
h.d <- c(253, 337, 395, 451, 495, 534, 573)
res.lm  <- lm(h.d ~ init.h)
res.lm2 <- update(res.lm,  . ~ . + I(init.h^2))
res.lm3 <- update(res.lm2, . ~ . + I(init.h^3))


## @knitr 
short_summary(res.lm2)


## @knitr 
alpha <- 0.05
tstar <- qt(1 - alpha/2, df=4)          # n=7; p=2; df=n-(p+1)
beta1 <- 1.05
SE <- 0.141
beta1 + c(-1,1 ) * tstar * SE


## @knitr 
anova(res.lm2,res.lm3)


## @knitr fig.keep="none"
pairs(stud.recs)


## @knitr 
d <- subset(stud.recs, select=-letter.grade)
res.lm <- lm(num.grade ~ ., data = d)
res.lm


## @knitr 
short_summary(res.lm)


## @knitr 
library(MASS)                 # load in MASS package for stepAIC
stepAIC(res.lm, trace=0)      # trace=0 suppresses intermediate output


## @knitr 
init.h <- c(600,700,800,950,1100,1300,1500)
h.d <- c(253, 337, 395, 451, 495, 534, 573)
res.lm3 <- lm(h.d ~ init.h + I(init.h^2) + I(init.h^3))
res.lm4 <- update(res.lm3, . ~ . + I(init.h^4))
anova(res.lm3, res.lm4)


## @knitr 
res <- lm(Volume ~ Girth + Height, data=trees)


## @knitr fig.keep="none"
plot(fitted(res),resid(res))      


## @knitr fig.keep="none"
res <- lm(attendance ~ year + runs.scored + wins + games.behind, 
          data=MLBattend)
plot(fitted(res), resid(res))


## @knitr 
short_summary(res)


## @knitr fig.keep="none"
res.1 <- lm(Deflection ~ Load, data=deflection)
res.2 <- update(res.1, . ~ . + I(Load^2))
plot(fitted(res.1), resid(res.1)) # clearly a bad model
plot(fitted(res.2), resid(res.2)) # looks good


## @knitr 
res.1 <- lm(weight ~ age + height, data=kid.weights)
res.2 <- update(res.1, . ~ . + I(height^2))
res.3 <- update(res.2, . ~ . + I(height^3))
res.4 <- update(res.3, . ~ . + I(height^4))
anova(res.1,res.2,res.3,res.4)


## @knitr 
res <- lm(body.fat ~ age + weight + height + BMI + neck +
          chest + abdomen + hip + thigh + knee + ankle + bicep + forearm + 
          wrist, data = fat)


## @knitr 
stepAIC(res, trace=0)            # after library(MASS)


## @knitr 
res <- lm(body.fat ~ ., data=fat[,-c(1,3,4)])      


## @knitr 
library(MASS)                 # load data set
res = lm(MPG.city ~ EngineSize + Weight + Passengers + Price, data=Cars93)
short_summary(res)
stepAIC(res, trace=0)


## @knitr 
x <- 1:10
y <- rnorm(10, 1 + 2*x + 3*x^2, 4)
require(MASS)
stepAIC(lm(y ~ x + I(x^2)), trace=0)       


## @knitr 
d <- with(baycheck, {
  n <- length(year)
  yt <- log(Nt[-1]/Nt[-n])
  nt <- Nt[-n]
  data.frame(yt, nt)
})


## @knitr 
d <- with(baycheck, {
  n <- length(year)
  yt <- log(Nt[-1]/Nt[-n])
  nt <- Nt[-n]
  data.frame(yt, nt)
})
lm(yt ~ nt, data=d)


