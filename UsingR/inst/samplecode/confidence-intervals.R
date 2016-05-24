
## @knitr universe-estiamtes, echo=FALSE, out.width=singlewide

  k <- 1:16
  labels <- age.universe$year
  rng <- age.universe[,1:2]
  
  plot.new()
  plot.window(xlim = range(k),ylim= range(rng,na.rm=T))
  axis(1,at=k,labels=labels, las=2)
  axis(2)
  title(main="Estimates for age of universe",
        xlab = "year of estimate",
        ylab = "Age (billions of years)")
  
  ## draw confidence intervals
  y.lower = rng$lower; y.upper = rng$upper
  
  lineDim = .2
  for(i in 1:length(k)) {
    if(!is.na(y.upper[i])) 
      lines(c(i-lineDim,i+lineDim),rep(y.upper[i],2))
    if(!is.na(y.lower[i]))
      lines(c(i-lineDim,i+lineDim),rep(y.lower[i],2))
  }
  
  y.lower[is.na(y.lower)] = min(rng,na.rm=T)
  y.upper[is.na(y.upper)] = max(rng,na.rm=T)
  
  
  
  for(i in 1:length(k)) {
    lines(c(i,i),c(y.lower[i],y.upper[i]),lty=8)
  }
  
  abline(h = 13.772, col=gray(.8))


## @knitr 
phat <- 0.46
n <- 1012
SE <- sqrt(phat*(1-phat)/n)
c(lower = phat - 1.96 * SE, upper = phat + 1.96 * SE)


## @knitr 
conf.level <- c(0.80, 0.90, 0.95, 0.99)
alpha <- 1 - conf.level
multipliers <- qnorm(1 - alpha/2)
setNames(multipliers, conf.level)       # give some names


## @knitr zstar-alpha, echo=FALSE, out.width=singlewide
## a picture showing z^* and alpha
the.cex = 1.5
## make a picture of density, and area between a and b

f = function(x) dnorm(x)

x1 = seq(-4,-1.5,length=500)
x2 = seq(-1.5,1.5,length=500)
x3 = seq(1.5,4,length=500)

## empty plot

plot(0,0,pch=" ",
     bty="n",
     xlim = c(-3,3),
     ylim = c(-.025,0.40),
     
     xlab = "x",
     ylab = "Normal density"
     )

shade.plot = function(x,y,col=jvgray) {
  ## x, y give a line, we want to make a polygon
  polygon(c(x[1],x,x[length(x)],x[1]),
          c(0,y,0,0),
          col=col)
}

shade.plot(x1,f(x1),col=gray(.7))
shade.plot(x2,f(x2),col=gray(1))
shade.plot(x3,f(x3),col=gray(.7))

## label 1-alpha
text(0,0.11,label=expression(1-alpha),cex=the.cex)
## label alpha/2's

text(-2.7,0.11,label=expression(alpha/2),cex=the.cex)
arrows(-2.5,0.09, -2.1, 0.05,angle = 45,length=.15)

text(2.6,0.11,label=expression(alpha/2),cex=the.cex)
arrows(2.5,0.09, 2.1, 0.05,angle = 45,length=.15)


## label z^*
text(-1.5,-.02,expression(-z^"*"),cex=the.cex)
text(1.5,-.02,expression(z^"*"),cex=the.cex)



## @knitr 
x <- 80; n <- 125
phat <- x/n
alpha <- 1 - 0.90
zstar <- qnorm(1 - alpha/2)
SE <- sqrt(phat * (1 - phat) / n)
MOE <- zstar * SE
phat + c(-1, 1) * MOE


## @knitr 
prop.test(x, n)


## @knitr 
binom.test(x, n)$conf.int


## @knitr 
confint(binom.test(x, n))


## @knitr 
zstar <- 0.03 / sqrt(.57*(1-.57)/1000)
alpha <- 2* pnorm(-zstar)
(1-alpha) * 100


## @knitr 
n <- 100; phat <- 0.45
prop.test(n*phat, n, conf.level = 0.8)
prop.test(n*phat, n, conf.level = 0.9)


## @knitr 
prop.test(5,100)


## @knitr 
prop.test(x=9, n=10, conf.level=0.80)$conf.int


## @knitr 
binom.test(x=9, n=10, conf.level=0.80)$conf.int


## @knitr 
phat <- .54
zstar <- 1.96
(zstar * sqrt(phat*(1-phat)) / 0.02)^2


## @knitr 
n <- 5;    SE <- sqrt(.8*(.2)/n); qnorm(.95) * SE
n <- 100;  SE <- sqrt(.8*(.2)/n); qnorm(.95) * SE
n <- 1000; SE <- sqrt(.8*(.2)/n); qnorm(.95) * SE


## @knitr 
n <- 250;  SE <- sqrt(phat*(1-phat)/n);  qnorm(.975) * SE
n <- 1000; SE <- sqrt(phat*(1-phat)/n); qnorm(.975) * SE


## @knitr 
## for 95% confidence
alpha <- 0.05; zstar <- qnorm(1 - alpha/2)
(zstar/0.01)^2/4
## for 90% confidence
alpha <- 0.1; zstar <- qnorm(1 - alpha/2)
(zstar/0.01)^2/4
## for 80% confidence
alpha <- 0.2; zstar <- qnorm(1 - alpha/2)
(zstar/0.01)^2/4


## @knitr 
M <- 50; n <- 20; p <- .5;     # toss 20 coins 50 times,
alpha <- 0.10;
zstar <- qnorm(1-alpha/2)
phat <- rbinom(M, n, p)/n        # divide by n for proportions
SE <- sqrt(phat*(1-phat)/n)    # compute SE


## @knitr 
sum(phat - zstar*SE < p & p < phat + zstar * SE)/n   


## @knitr eval=FALSE
## matplot(rbind(phat - zstar*SE, phat + zstar*SE),
##          rbind(1:m,1:m),type="l",lty=1)
## abline(v=p)                  # indicate parameter value


## @knitr echo=FALSE
alpha <- 0.95


## @knitr 
tstar <- qt(1 - alpha / 2, df=n - 1)
alpha <- 2 * pt(-tstar, df=n - 1)    


## @knitr 
zstar <- qnorm(1 - alpha / 2)
alpha <- 2 * pnorm(-zstar)    


## @knitr 
xbar <- 66; s <- 4; n <- 30
alpha <- 1 - 0.80
tstar <- qt(1 - alpha/2, df = n-1)      # 1.311
SE <- s/sqrt(n)
MOE <- tstar * SE                      
xbar + c(-1,1) * MOE


## @knitr fig.keep="none"
ozs <- c(1.95, 1.80, 2.10, 1.82, 1.75, 2.01, 1.83, 1.90)
qqnorm(ozs)                             # approximately linear
confint(t.test(ozs, conf.level=0.80))


## @knitr echo=FALSE
lower <- round(t.test(ozs, conf.level=0.80)$conf.int[1], 2)
upper <- round(t.test(ozs, conf.level=0.80)$conf.int[2], 2)


## @knitr echo=FALSE
set.seed(10)


## @knitr 
k <- 5                                  # answers per question
n <- 10                                 # number of questions
tstat <- function(x, mu=0) (mean(x) - mu)/(sd(x)/sqrt(length(x)))
res <- replicate(2000, {
  xs <- sample(1:k, n, replace=TRUE)
  tstat(xs, (k+1)/2)
})
## look at one quantile
sum(res > qt(0.05, df=n-1, lower.tail=FALSE)) / length(res)


## @knitr echo=FALSE
n <- 10
tstat <- function(x, mu=0) (mean(x) - mu)/(sd(x)/sqrt(length(x)))
l <- lapply(c(2,5,15), function(k) replicate(2000, {tstat(sample(1:k, n, replace=TRUE), (k+1)/2)}))


## @knitr echo=FALSE, out.width=triplewide
ts <- rt(2000, df=n-1)
qqplot(l[[1]], ts, main="k=2")
qqplot(l[[2]], ts, main="k=5")
qqplot(l[[3]], ts, main="k=15")


## @knitr 
x <- c(175, 185, 170, 184, 175)
t.test(x,conf.level = 0.90, alt="less")


## @knitr 
xbar <- 9.5; s <- 1; n <- 15
alpha <- 0.1
tstar <- qt(1 - alpha/2, df=n-1)
SE <- s/sqrt(n)
c(xbar - tstar*SE, xbar + tstar*SE)


## @knitr 
hist(stud.recs$sat.m)
t.test(stud.recs$sat.m, conf.level=0.80)


## @knitr 
t.test(homedata$y1970, conf.level = 0.9)
t.test(homedata$y2000, conf.level = 0.9)


## @knitr eval=FALSE
## qqnorm(homedata$y2000)                  # no, skewed
## length(homedata$y2000)                  # yes, n=6841


## @knitr 
yr5 <- subset(kid.weights, subset= 5*12 <= age & age < 6*12)


## @knitr fig.keep="none"
hist(yr5$weight)                        # not long tailed, but skewed
length(yr5$weight)                      # not large -- may be issues!
t.test(yr5$weight, conf.level=0.90)


## @knitr 
hist(brightness)              # bell-shaped, looks good
t.test(brightness, conf.level = 0.90)$conf.int


## @knitr 
t.test(normtemp$temperature,conf.level=.90)


## @knitr 
finger <- with(Macdonell, rep(finger, frequency))


## @knitr echo=FALSE
set.seed(10)


## @knitr 
require(HistData)
finger <- with(Macdonell, rep(finger, frequency))
res <- replicate(750, mean(sample(finger, 4, replace=TRUE)))
quantile(res, c(0.025, 0.975))


## @knitr 
n <- 4
x <- sample(finger, n, replace=TRUE)
t.test(x, conf.level=0.95)$conf.int


## @knitr fig.keep="none"
n <- 10; M <- 2000
res <- replicate(M, {
  x <- rnorm(n)
  (mean(x) - 0)/(sd(x)/sqrt(n))
})
qqplot(res, rt(M, df=n-1))              # compare to t-distribution


## @knitr fig.keep="none"
n <- 10
## simulated data
boxplot(rt(1000, df=n - 1), rnorm(1000))
## theoretical qqplots
x <- seq(0, 1, length=150)
plot(qt(x, df=n - 1), qnorm(x))
abline(0, 1)
## compare densities
curve(dnorm(x), -3.5, 3.5)
curve(dt(x, df=n - 1), lty=2, add=TRUE)


## @knitr ech0=FALSE
set.seed(10)


## @knitr 
M <- 200; n <- 10
alpha <- 1 - 0.90
res <- replicate(M, {
  x <- rnorm(n, mean=0, sd=2)
  c(z=qnorm(1-alpha/2) * 2/sqrt(n), 
    t=qt(1 - alpha/2, df=n-1) * sd(x)/sqrt(n))
})                                      # a matrix
sum(res["z",] < res["t", ]) / M


## @knitr 
n <- 10; alpha <- 1 - 0.90              # say
lstar <- qchisq(alpha/2, df=n-1)
rstar <- qchisq(1-alpha/2, df=n-1)


## @knitr 
s2 <- 12; n <- 10
alpha <- 1 - 0.95
lstar = qchisq(alpha/2, df=n - 1)
rstar = qchisq(1 - alpha/2, df=n - 1)
(n-1) * s2 * c(1/rstar, 1/lstar)        # CI for sigma squared


## @knitr 
n <- 11; m <- 16
alpha <- 1 - 0.90
qf(c(alpha/2, 1- alpha/2), df1=n-1, df2=m-1)


## @knitr 
n <- 10; m <- 20
alpha <- 1 - 0.80
lr <- qf(c(alpha/2, 1- alpha/2), df1=n-1, df2=m-1)
lr
s <- (2.3/2.8)^2
sqrt(s/rev(lr))


## @knitr 
place <- nym.2002$place
n <- length(place)
alpha <- 1 - 0.90
(1-alpha)^(1/n)
max(place)/(1-alpha)^(1/n)
max(place)


## @knitr obama-approval, echo=FALSE,out.height=".75\\textwidth", out.width=".98\\textwidth"
x <- ObamaApproval
x$start <- as.Date(x$start)
orgs <- names(rev(sort(table(x$org))))

plot(approve ~ start, x, xlab="date", col=NA, bty="L")


cols <- gray(seq(.25, .75, length=length(orgs)), alpha=.75)
names(cols) <- orgs



for(i in 1:nrow(x)) {
  p <- x$approve[i]/100
  d <- x$start[i]
  n <- x$n[i]
  points(d, 100*p, col=cols[x$org[i]], pch=16)

  moe = 1.96 * sqrt(p*(1-p)/n)
  lines(c(d,d),
        100*p + c(-moe, moe) * 100,
        lty=2,
        col=cols[x$org[i]])
}



## @knitr 
prop.test(x=c(560,570), n=c(1000,1200), conf.level=0.95)


## @knitr drug-placebo-boxplot, echo=FALSE, fig.keep="none"
x <- c(0, 0, 0, 2, 4, 5, 13, 14, 14, 14, 15, 17, 17)    
y <- c(0, 6, 7, 8, 11, 13, 16, 16, 16, 17, 18)
boxplot(list(placebo=x, ephedra=y), col="gray") # compare spreads


## @knitr echo=FALSE, out.width=doublewide
x <- c(0, 0, 0, 2, 4, 5, 13, 14, 14, 14, 15, 17, 17)    
y <- c(0, 6, 7, 8, 11, 13, 16, 16, 16, 17, 18)
boxplot(list(placebo=x, ephedra=y), col="gray") # compare spreads


## @knitr fig.keep="none"
x <- c(0, 0, 0, 2, 4, 5, 13, 14, 14, 14, 15, 17, 17)    
y <- c(0, 6, 7, 8, 11, 13, 16, 16, 16, 17, 18)
boxplot(list(placebo=x, ephedra=y), col="gray") # compare spreads


## @knitr 
confint(t.test(x,y, var.equal=TRUE))


## @knitr 
confint(t.test(x,y))


## @knitr 
out <- t.test(childHeight ~ gender, GaltonFamilies, conf.level=0.95)
confint(out)


## @knitr 
library(MASS)
names(shoes)
### Alternately: with(shoes, t.test(A,B,conf.level=0.9,paired=TRUE))
with(shoes, {
  out <- t.test(A-B, conf.level=0.9)
  confint(out)
})


## @knitr 
c1 <- c(3.1, 3.3, 1.7, 1.2, 0.7, 2.3, 2.9)
c2 <- c(1.8, 2.3, 2.2, 3.5, 1.7, 1.6, 1.4)
t.test(c1,c2, conf.level=0.80)


## @knitr 
t.test(c1,c2, conf.level=0.80, var.equal=TRUE)$conf.int


## @knitr fig.keep="none"
gp.400 <- c(7,0,8,1,10,12,2,9,5,2)
gp.1200 <- c(2,1,5,1,4,7,-1,8,7,3)
boxplot(list(gp.400=gp.400, gp.1200 = gp.1200))


## @knitr 
t.test(gp.400, gp.1200, conf.level = 0.90, var.equal=FALSE)


## @knitr fig.keep="none"
foster <- c(80, 88,75, 113, 95, 82, 97, 94, 132, 108)
biological <- c(90, 91, 79, 97, 97, 82, 87, 94, 131, 115)
plot(foster, biological)
boxplot(foster - biological)
t.test(foster, biological, conf.level=0.90, paired=TRUE)      


## @knitr 
plot(age ~ dage, data=babies)
b <- subset(babies, subset= dage < 99 & age < 99 ) # fix 99 = NA
plot(age ~ dage, data=b)
hist(b$dage - b$age)
t.test(b$dage,b$age, conf.level=0.95, paired=TRUE)


## @knitr echo=FALSE
set.seed(10)


## @knitr 
ceo <- ceo2013[sample(1:200, 10),
               c("company_name_ticker", "cash_compensation")]


## @knitr 
ceos <- round(ceo[,2]/1e6, 2)           # in millions
n <- length(ceos)
pbinom(0:n,n,1/2)         


## @knitr 
sort(ceos)


## @knitr 
alpha <- 1 - 0.90
n <- length(ceos)
j <- qbinom(alpha/2, n, 1/2)
sort(ceos)[c(j, n+1-j)] 


## @knitr 
ans <- wilcox.test(log(ceos), conf.int=TRUE, conf.level=0.9)
confint(ans)
confint(ans, transform=exp)             # inverse of log.


## @knitr 
ceos_past <- ceo2013[sample(1:200, 10), "cash_compensation_past"]
ceos_past <- round(ceos_past / 1e6, 2)  # in millions
ceos_past
ceos                                     # 2013 data


## @knitr ceo-trough, echo=FALSE, out.width=singlewide
d1 <- density(ceos); d2 <- density(ceos_past)
xrng <- range(c(d1$x, d2$x))
yrng <- range(c(d1$y, d2$y))
plot(d1, xlim=xrng, ylim=yrng, lty=1, main="CEO compensation")
lines(d2, lty=2)
legend(.25*xrng[1] + .75*xrng[2], .25*yrng[1] + .75*yrng[2], 
       legend=c("CEO 2013", "CEO past"), lty=1:2)


## @knitr 
wilcox.test(ceos, ceos_past, conf.int=TRUE, conf.level=0.9)


## @knitr 
commutes <- rep(c(21,22,23,24,25,26,29,31,33),
  c(3,2,2,6,2,1,1,2,1))
ans <- wilcox.test(log(commutes), conf.int=TRUE, conf.level=0.9)
confint(ans)
confint(ans, transform=exp)             # also exp(ans$conf.int)


## @knitr 
wilcox.test(u2$October,u2$"The Joshua Tree", conf.int=TRUE)


## @knitr fig.keep="none"
hist(cfb$AGE)                 # symmetric
wilcox.test(cfb$AGE, conf.int=TRUE)


## @knitr fig.keep="none"
hist(log(cfb$INCOME +1))
ans <- wilcox.test((log(cfb$INCOME +1)), conf.int=TRUE)
exp(ans$conf.int) - 1


## @knitr 
n <- 20
M <- 250
res <- replicate(M, {
  x <- rnorm(n)
  sum(rank(abs(x))[x>0])                # only add positive values
})


## @knitr fig.keep="none"
hist(res, probability=TRUE)
x <- 40:140
lines(x, dsignrank(x,n))                # density-like, but discrete


## @knitr fig.keep="none"
qqplot(res,rsignrank(100,n))    


