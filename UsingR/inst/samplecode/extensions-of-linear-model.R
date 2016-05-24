## display a short summary of a linear model
short_summary <- function(x) {
  x <- summary(x)
  cmat <- coef(x)
  printCoefmat(cmat)
}


## @knitr echo=FALSE, out.width=singlewide
logisticfn <- function(x) exp(x) / (1 + exp(x))
curve(logisticfn, -5, 5, main="Logistic function")
points(0, logisticfn(0), pch=15, cex=2)


## @knitr 
x1 <- rep(1:10, 2)
x2 <- rchisq(20, df=2)
y <- rnorm(20, mean=x1 + 2*x2, sd=2)
res.lm <- lm(y ~ x1 + x2)
short_summary(res.lm)


## @knitr 
res.glm <- glm(y ~ x1 + x2, family="gaussian")
summary(res.glm)


## @knitr 
babies.prem = subset(babies,
  subset= gestation < 999 & wt1 < 999 & ht < 99 & smoke < 9, 
  select=c("gestation","smoke","wt1","ht"))    


## @knitr 
babies.prem$preemie = with(babies.prem, as.numeric(gestation < 7*37))
table(babies.prem$preemie)


## @knitr fig.keep="none"
babies.prem$BMI = with(babies.prem, (wt1/2.2) / (ht*2.54/100)^2)
hist(babies.prem$BMI)                     # looks okay


## @knitr 
res <- glm(preemie ~ factor(smoke) + BMI, family=binomial,
           data=babies.prem)
short_summary(res)


## @knitr 
library(MASS)
stepAIC(res, trace=0)


## @knitr 
first.name <- gl(2, 2500, 5000, labels=c("yes", "no"))
offer      <- gl(2, 1250, 5000, labels=c("yes", "no"))
opened <- c(rep(1:0, c(20, 1250-20)), rep(1:0, c(15, 1250-15)),
            rep(1:0, c(17, 1250-17)), rep(1:0, c( 8, 1250-8)))
xtabs(opened ~ first.name + offer)


## @knitr 
f <- function(x) rep(1:0, c(x, 1250-x))  
opened <- c(sapply(c(20, 15, 17, 8), f))


## @knitr 
res.glm <- glm(opened ~ first.name + offer, family="binomial")
short_summary(res.glm)


## @knitr 
opened <- c(8,15,17,20)
opened.mat <- cbind(opened=opened, not.opened=1250 - opened)
opened.mat


## @knitr 
offer <- c(0, 0, 1, 1)
first.name <- c(0, 1, 0, 1)


## @knitr 
glm(opened.mat ~ first.name + offer, family="binomial")


## @knitr yellow-fin-plot, fig.keep="none"
plot(count ~ year, data=yellowfin)


## @knitr 
f <- function(t, N, r, d) N*(exp(-r*(t-1952))*(1-d) + d)  


## @knitr yellow-fin-n-6, fig.keep="none", eval=FALSE
## curve(f(x, N=6, r=1/10, d=0.1), add=TRUE)


## @knitr 
res.yf <- nls(count ~ f(year, N, r, d), start=c(N=6, r=1/10, d=0.1), 
  data=yellowfin)
res.yf


## @knitr add-exponential, fig.keep="none", eval=FALSE
## curve(f(x, N=6.02, r=0.0939, d=0.0539), add=TRUE, lty=2, lwd=2)
## legend(1980, 6, legend=c("exploratory", "exponential"), lty=1:2)


## @knitr echo=FALSE, out.width=singlewide
plot(count ~ year, data=yellowfin)
  title(main="Mean catch per 100 hooks")
curve(f(x, N=6, r=1/10, d=0.1), add=TRUE)  
curve(f(x, N=6.02, r=0.0939, d=0.0539), add=TRUE, lty=2, lwd=2)  
legend(1980, 6, legend=c("exploratory", "exponential"), lty=1:2)


## @knitr eval=FALSE
## tmp <- 1952:2000
## lines(tmp, predict(res.yf, data.frame(year = tmp)))


## @knitr urchin-growth, fig.keep="none"
logistic <- function(t, Y, k, t0, m) Y * (1 + exp(-k * (t-t0)))^(-1)
richards <- function(t, Y, k, t0, m) Y * (1 - exp(-k * (t-t0)))^m

plot(jitter(size) ~ jitter(age,3), data=urchin.growth,
     xlab="age", ylab="size", main="Urchin growth by age")


## @knitr fig.keep="none", eval=FALSE
## curve(logistic(x, Y=60, k=1, t0=2), add=TRUE)


## @knitr 
res.logistic <- nls(size ~ logistic(age, Y, k, t0), 
                    start = c(Y=60, k=1, t0=2),
                    data  = urchin.growth)
res.logistic


## @knitr eval=FALSE
## curve(logistic(x, Y=53.903, k=1.393, t0=1.958), add=TRUE)


## @knitr 
AIC(res.logistic)


## @knitr eval=FALSE
## curve(richards(x, Y=53.903, k=1.393, t0=1.958, m = 1), add=TRUE)
## legend(4, 20, legend=c("logistic growth", "Richards"), lty=1:2)


## @knitr 
res.logistic <- nls(size ~ logistic(age,Y,k,t0,m), 
                    data  = urchin.growth,
                    start = c(Y=53, k=1.393, t0=1.958, m=1))


## @knitr richardson, fig.keep="none"
res.richards <- nls(size ~ richards(age, Y, k, t0, m), 
                    data  = urchin.growth,
                    start = c(Y=53, k=0.5, t0=0, m=1))
res.richards


## @knitr eval=FALSE
## curve(richards(x, Y=57.26, k=0.78, t0=-0.8587, m = 6.0636),
##       add=TRUE, lty=2)


## @knitr 
AIC(res.richards)


## @knitr echo=FALSE
plot(jitter(size) ~ jitter(age,3), data=urchin.growth,
     xlab="age", ylab="size", main="Urchin growth by age")
curve(logistic(x, Y=53.903, k=1.393, t0=1.958),             add=TRUE)    
curve(richards(x, Y=57.26, k=0.78, t0=-0.8587, m = 6.0636), add=TRUE, lty=2)
legend(5, 30, lty=1:2, legend=c("Logistic", "Richards"))


## @knitr 
res <- glm(enjoyed ~ gender + age, data=tastesgreat, 
           family=binomial)
summary(res)


## @knitr 
library(MASS)                 # load stepAIC
stepAIC(glm(healthy ~ p + g, healthy, family=binomial))


## @knitr 
library(MASS)
res.glm <- glm(low ~ age + lwt + smoke + ht + ui, data=birthwt,
 family = binomial)
summary(res.glm)


## @knitr 
stepAIC(res.glm, trace=0)


## @knitr 
hfm = hall.fame$Hall.Fame.Membership != "not a member"    


## @knitr 
hfm <- hall.fame$Hall.Fame.Membership != "not a member"    
res <- glm(hfm ~ BA + HR + hits + games, data=hall.fame, 
  family="binomial")
stepAIC(res, trace=0)


## @knitr 
res.full <- glm(cbind(ncases, ncontrols) ~ agegp + tobgp * alcgp,
                data = esoph, family = binomial())


## @knitr 
res.add  <- glm(cbind(ncases, ncontrols) ~ agegp + tobgp + alcgp,
                 data = esoph, family = binomial())


## @knitr 
res.full <- glm(cbind(ncases, ncontrols) ~ agegp + tobgp * alcgp,
 data = esoph, family = binomial())
res.add  <- glm(cbind(ncases, ncontrols) ~ agegp + tobgp +  alcgp,
 data = esoph, family = binomial())
AIC(res.full)
AIC(res.add)


## @knitr fig.keep="none"
plot(circumference ~ age, data=Orange, subset = Tree == 1)


## @knitr 
g <- function(t, Y, k, t_0) Y*( 1 + exp(-k*(t-t_0)))^(-1)
res.tree <- nls(circumference ~ g(age, Y, k, t0),
                start=c(Y=150, k=1/300, t0=750), 
                data=Orange, subset=Tree == 1)
res.tree


## @knitr cache=FALSE
age <- 0:1500
plot(circumference ~ age, data=Orange, subset = Tree == 1)
lines(age, predict(res.tree, data.frame(age=age)))      


## @knitr fig.keep="none"
f <- function(t,Y,k,t0)  Y * (1 + exp(-k*(t-t0)))^(-1)
plot(weight ~ Time, data=ChickWeight, subset= Chick == 1) 


## @knitr 
nls(weight ~ f(Time, Y, k, t0), start=c(Y=400, k=1/10, t0=20), 
 subset= Chick == 1, data=ChickWeight)


## @knitr fig.keep="none"
curve(f(x, Y=937, k=.08768, t0=35.22270), add=T)      


## @knitr fig.keep="none"
library(MASS)                 # load in help page
example(wtloss) 
wtloss.fm
plot(Weight ~ Days, wtloss)   # make scatterplot
lines(predict(wtloss.fm) ~ Days, data = wtloss)


## @knitr 
predict(wtloss.fm, newdata=data.frame(Days=365))


## @knitr 
l  <- function(t, a, b, k, t0) (a + b*t) * (1 - exp(-k*(t-t0)))
l1 <- function(t, a, k, t0) l(t, a, 0, k, t0)
res.l <- nls(length ~ l(age,a,b,k,t0), data=reddrum,
             start=c(a=32,b=.25,k=.5,t0=0))
res.l1 <- nls(length ~ l1(age,a,k,t0), data=reddrum,
              start <- c(a=32,k=.5,t0=0))
AIC(res.l)
AIC(res.l1) 


## @knitr 
year <- with(midsize, 2004-Year)   
f <- function(x, N, r) N * exp(-r*x)
with(midsize, nls(Accord ~ f(year,N,r), start=c(N=Accord[1], r=1/5)))
with(midsize, nls(Camry ~ f(year,N,r),  start=c(N=Camry[1],  r=1/5)))
with(midsize, nls(Taurus ~ f(year,N,r), start=c(N=Taurus[1], r=1/5)))


