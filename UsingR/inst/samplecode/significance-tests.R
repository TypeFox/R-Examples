
## @knitr echo=FALSE, out.width=doublewide
par(cex.main=1.5)
par(cex.lab=1.5)
par(cex.sub=1.5)

shadenorm <- function (obs = 0.7, std=1) {
xmin <- -2; xmax<-3

ho <- 0
ha <- 1

colo <- rgb(.5, .5, .5, .5)
cola <- rgb(.7, .7, .7, .5)

sa<-so<-std


x <- seq(xmin, xmax, length=200)

plot.new()
plot.window(xlim=c(xmin, xmax), ylim=c(0, dnorm(0, sd=sa)))

lines(x, dnorm(x, mean=ho, sd=so), lwd=2)
lines(x, dnorm(x, mean=ha, sd=sa), col=cola, lwd=2)

x <- seq(obs, xmax, length=200)
polygon(c(x, rev(x)), c(dnorm(x, mean=ho, sd=so), 0*x), col=colo)

x <- seq(xmin, obs, length=200)
polygon(c(x, rev(x)), c(dnorm(x, mean=ha, sd=sa), 0*x), col=cola)

axis(1, at=c(-2, -1, 0, 0.7, 1, 2, 3))

}

shadenorm()
shadenorm(std=1/sqrt(10))



## @knitr echo=FALSE, out.width=triplewide
par(cex.main=1.5)
par(cex.lab=1.5)
par(cex.sub=1.5)
shadeplot <- function(obs=1.5, left=TRUE, right=TRUE, main="") {
  plot.new()
  plot.window(xlim=c(-3,3), ylim=c(0, dnorm(0)))
  
  col=rgb(.7, .7, .7, .5)
  
  x <- seq(-3, 3, length=200)
  lines(x, dnorm(x, 0, 1), lwd=2)
  
  if(left) {
    x <- seq(-3, -abs(obs), length=200)
    polygon(c(x, rev(x)), c(dnorm(x, 0, 1), 0*x), col=col)
  }
  
  if(right) {
    x <- seq(abs(obs), 3, length=200)
    polygon(c(x, rev(x)), c(dnorm(x, 0, 1), 0*x), col=col)
  }
  abline(h=0)
  
  title(main=main)
  
}
  
obs = 1.5
shadeplot(obs, left=TRUE, right=FALSE, main="less")
text(-obs, 0, "observed", pos=1)

shadeplot(obs, left=FALSE, right=TRUE, main="greater")
text(obs, 0, "observed", pos=1)

shadeplot(obs, left=TRUE, right=TRUE, main="two-sided")
text(obs, 0, "observed", pos=1)



## @knitr 
phat <- 22695 / 150000                   # 0.1513
p0 <- 0.1500; n <- 150000
SD <- sqrt(p0 * (1-p0)/n)
pnorm(phat, mean=p0, sd=SD, lower.tail=FALSE) # p-value


## @knitr 
prop.test(x=22695, n=150000, p=.1500, alt="greater")


## @knitr 
out <- prop.test(x=22695, n=150000, p=.1500, alt="two.sided")
out$p.value                             # just the p-value


## @knitr 
table(samhda$marijuana)


## @knitr 
x <- 309; n <- 309 + 281
prop.test(x, n, p=.5, alt="greater")


## @knitr 
out <- prop.test(40,50, p = .75, alt="greater")
out$p.value


## @knitr 
out <- prop.test(40, 75, p = 0.281, alt="two.sided")
out$p.value


## @knitr 
out <- prop.test(5731, 5760, p=0.999, alt="less")
out$p.value


## @knitr 
p <- 0.1500; n <- 150000
out <- qnorm(.95, mean=p,sd=sqrt(p*(1-p)/n))


## @knitr 
n * out


## @knitr 
out <- prop.test(x=2700, n=25000, p=0.1, alt="greater" )
out$p.value


## @knitr 
out <- prop.test(x=.16*200, n=200, alt="two.sided")
out$p.value


## @knitr echo=FALSE
SUV <- c(11.4, 13.1, 14.7, 14.7, 15.0, 15.5, 15.6, 15.9, 16.0, 16.8)
stem(SUV, scale=2)


## @knitr 
mpg <- c(11.4, 13.1, 14.7, 14.7, 15.0, 15.5, 15.6, 15.9, 16.0, 16.8)
xbar <- mean(mpg); s <- sd(mpg); n <- length(mpg)
c(xbar=xbar, s=s, n-n)                  # summaries
SE <- s/sqrt(n)
obs <- (xbar - 17)/SE
pt(obs, df = 9, lower.tail = T)


## @knitr 
t.test(mpg, mu = 17, alt="less")


## @knitr 
costs <- c(304, 431, 385, 987, 303, 480, 455, 724, 642, 506)


## @knitr 
stem(costs)


## @knitr 
t.test(costs, mu=500, alt="greater")


## @knitr echo=FALSE, out.width=triplewide
par(cex.main=1.5)
par(cex.lab=1.5)
par(cex.sub=1.5)

show_power <- function (n = 1, std=1) {
  xmin <- -2; xmax<-3

  ho <- 0
  ha <- 1
  
  colo <- rgb(.5, .5, .5, .5)
  cola <- rgb(.7, .7, .7, .5)
  
  sa<-so<-std/sqrt(n)
  
  alpha <- 0.05
  obs <- qnorm(1 - alpha, mean=ho, sd=so)
  
  x <- seq(xmin, xmax, length=200)
  
  plot.new()
  plot.window(xlim=c(xmin, xmax), ylim=c(0, dnorm(0, sd=sa)))
  title(main=sprintf("P(Type-II error) = %.2f", pnorm(obs, mean=ha, sd=sa)),
        xlab=sprintf("n=%s", n))
  
  
  lines(x, dnorm(x, mean=ho, sd=so), lwd=2)
  lines(x, dnorm(x, mean=ha, sd=sa), col=cola, lwd=2)
  
  x <- seq(obs, xmax, length=200)
  polygon(c(x, rev(x)), c(dnorm(x, mean=ho, sd=so), 0*x), col=colo)
  
  x <- seq(xmin, obs, length=200)
  polygon(c(x, rev(x)), c(dnorm(x, mean=ha, sd=sa), 0*x), col=cola)
  
  axis(1, at=c(-2, -1, 0, 0.7, 1, 2, 3))
  
}
show_power(n=1)
show_power(n=5)
show_power(n=10)


## @knitr power.t.test
alpha <- 0.05; beta <- 0.20
power.t.test(delta=1, sd=1, 
             sig.level=alpha, power=1-beta, 
             type="one.sample", alternative="one.sided")


## @knitr 
xbar = 58260; n = 25; sd = 3250
mu = 55000; SE = sd/sqrt(n)
T = (xbar - mu)/SE
T
pt(T, df=n-1, lower.tail=FALSE)


## @knitr 
xbar <- 4.03; s <- 0.42; n <- 800
mu <- 4.00
SE <- s/sqrt(n)
obs <- (xbar - mu)/SE
pnorm(obs, lower.tail=FALSE)


## @knitr fig.keep="none"
hist(stud.recs$sat.m)         # n is large, data is short-tailed
out <- t.test(stud.recs$sat.m, mu=500, alt="two.sided")
out$p.value <= 0.05                     # TRUE, reject


## @knitr 
x <- babies$dht
x <- x[x<99]
t.test(x, mu=68, alt="greater")


## @knitr 
xbar <- 0.5; s <- 3.77; n <- 7; 
mu <-  0
SE <- s/sqrt(n)
obs <- (xbar - mu)/SE
2*pt(obs, df = n-1, lower.tail=FALSE)


## @knitr fig.keep="none"
hist(OBP)
length(OBP)


## @knitr 
t.test(OBP, mu = 0.330, alt="two.sided")


## @knitr fig.keep="none"
x <- normtemp$temperature
hist(x)
t.test(x, mu = 98.6, alt="two.sided")


## @knitr 
alpha <- 0.05; beta <- 0.20
delta <- 1; std <- 1
power.t.test(sd=std, sig.level=alpha, power=1-beta, delta=delta,
             type="one.sample", alternative="one.sided")


## @knitr echo=FALSE, out.width=singlewide

jvgray=gray(.8)
## show duality picture
op = par(no.readonly=T)
par(mai=c(0,0.5698,0, 0.2919))
plot.norm.shade = function(x,mean=0,sd=1,col="transparent") {
  y = dnorm(x,mean,sd)
  polygon(c(x,rev(x)),c(y,rep(0,length(y))),col=col)
}

plot(0,0,pch= " ",
     xlim=c(-3,3),ylim=c(-.03,.42),
     xlab="",ylab="",
     xaxt="n",yaxt="n",
     bty="n"
     )


x1 = seq(-3,-1.96,length=100)
x2 = seq(-1.96,1.96,length=100)
x3 = seq(1.96,3,length=100)

plot.norm.shade(x1,col=jvgray)
lines(x2,dnorm(x2))
plot.norm.shade(x3,col=jvgray)
abline(h=0)

xbar = 2.5
SE = 1.96
under = -0.03
over  = 0.03
text(0,under,expression(mu),cex=1.5)
text(SE,under,expression(mu+t^"*" *SE),cex=1.5)
text(xbar,over*1.2,expression(bar(x)),cex=1.5)
text(xbar-SE,over*1.2,expression(bar(x)-t^"*" *SE),cex=1.5)
arrows(xbar-SE,over/2,xbar+SE,over/2,code=3,angle=0)
text(xbar-SE,over/2,"(")
points(xbar,over/2,pch=16)

text(0,.11,"rejection region")
arrows(-.5,.09,-1.9,.04,length=.15)
arrows(.5,.09,1.9,.04,length=.15)


par(op)



## @knitr 
calls <- c(2, 1, 3, 3, 3, 3, 1, 3, 16, 2, 2, 12, 20, 3, 1)
obs <- sum(calls > 5)          # find observed value of T
n <- length(calls)
1 - pbinom(n-obs - 1, n, 1/2)  # we want P(T >= 12) = 1 - P(T <= 11)


## @knitr 
k <- max(obs, n - obs)                  # k is 12
2*(1 - pbinom(k-1 , n, 1/2))


## @knitr 
wilcox.test(log(salmon.rate), mu=log(0.005), alt="greater")


## @knitr 
T <- sum(salmon.rate > .005); n <- length(salmon.rate)
1 - pbinom(T - 1, n, 1/2)


## @knitr 
T <- sum(exec.pay > 22)
n <- length(exec.pay)
T
pbinom(T - 1, n, .5, lower.tail=FALSE)


## @knitr 
wilcox.test(log(exec.pay), mu = log(22), alt="greater")


## @knitr fig.keep="none"
ep <- exec.pay[exec.pay > 0]
hist(log(ep))  


## @knitr 
smokers <- subset(babies, smoke == 1 & gestation != 999)


## @knitr 
smokers <- subset(babies, smoke == 1 & gestation != 999)
wilcox.test(smokers$gestation, mu = 40*7)


## @knitr echo=FALSE
set.seed(10)


## @knitr 
m <- 200; n <- 10

out <- replicate(m, {
  x <- rnorm(n, mean=1, sd=2)
  ttest <- t.test(x, mu=0, alt = "greater")$p.value
  sgntest <- 1 - pbinom(sum(x > 0) - 1, n, 1/2)
  c(t.test   = ifelse(ttest   < 0.05, 1, 0),
    sign.test= ifelse(sgntest < 0.05, 1, 0))
})

res.t <- out["t.test",]
res.sign <- out["sign.test",]
results <- c(t = sum(res.t)/m, sign=sum(res.sign)/m)


## @knitr 
phat <- c(0.1500, 0.1513)               # the sample proportions
n <- c(160000, 150000)                  # the sample sizes
n*phat                                  # the counts
prop.test(n*phat, n, alt="less")


## @knitr 
p <- sum(n*phat)/sum(n)        # (n_1p_1 + n_2p_2)/(n_1 + n_2)
obs <- (phat[1]-phat[2])/sqrt(p*(1-p)*sum(1/n))
obs 
pnorm(obs)


## @knitr 
x <- c("A"=14, "B"=15)
n <- c("A"=150, "B"=125)
prop.test(x, n, alt="less")


## @knitr 
n <- c(600, 1050)
x <- c(250, n[2] * 0.52)
prop.test(x, n)                         # use default alternative


## @knitr radio-station
phat <- c(0.82, 0.70)
n <- c(350, 350)
x <- n * phat
prop.test(x, n)                         # using default alternative


## @knitr 
x <- c("had"=153, "didnot"=196)
n <- c("had"=30000, "didnot"=30000)
prop.test(x, n)                         # use two.sided default


## @knitr 
x <- c(18, 3)
n <- c(22, 22)
prop.test(x, n)


## @knitr 
1250*(1-.99)


## @knitr 
percents <- c(98.9,96.9)
n <- c(1250, 1100)
p <- percents/100
x <- n * p
prop.test(x, n)                         # use two.sided default


## @knitr 
p <- c(0.31, 0.41)
n <- c(1000, 100)
x <- n * p
prop.test(x, n, alt="less")


## @knitr 
p <- c(.216,.193)
n <- c(10000,10000)
x <- n * p
out <- prop.test(x,n)
out$p.value < 0.01                      # TRUE


## @knitr 
m300 <- c(284, 279, 289, 292, 287, 295, 285, 279, 306, 298)  
m600 <- c(298, 307, 297, 279, 291, 335, 299, 300, 306, 291)


## @knitr m300-m600,eval=FALSE
## plot(density(m300))
## lines(density(m600), lty=2)


## @knitr 
t.test(m300, m600, var.equal=TRUE)


## @knitr 
t.test(m300, m600)


## @knitr echo=FALSE, out.width=doublewide
plot(density(m300))
lines(density(m600), lty=2)


## @knitr 
Finasteride <- c(5, 3, 5, 6, 4, 4, 7, 4, 3)
placebo <- c(2, 3, 2, 4, 2, 2, 3, 4, 2)
t.test(Finasteride, placebo, paired=TRUE, alt="two.sided")


## @knitr fig.keep="none"
pre  <- c(77, 56, 64, 60, 57, 53, 72, 62, 65, 66)
post <- c(88, 74, 83, 68, 58, 50, 67, 64, 74, 60)
boxplot(pre, post)
out <- t.test(pre, post, var.equal=TRUE, alt="less")
out$p.value


## @knitr 
out <- t.test(pre, post, paired=TRUE, alt="less")
out$p.value


## @knitr echo=FALSE, out.width=singlewide

## two sample figure
colo <- rgb(0, 0, 0, .5)
cola <- rgb(.4, .4, .4, .5)


xlim=c(0,20)

xplus = 4
df = 5
x = rchisq(8,df)
y = rchisq(8,df) + xplus


pts = seq(xlim[1],xlim[2],length=250)
c1 = dchisq(pts,5)
plot(pts,c1,type="l",xlim=xlim+c(0,xplus),
     bty="n",
     xlab = "",
     ylab=""
     )
points(pts+xplus,c1,lty=2,type="l")

points(x, 0*x, pch=16, col=colo)
points(y, 0*y + 0.01, pch=16, col=cola)



## @knitr fig.keep="none"
A <- c(5.8, 1.0, 1.1, 2.1, 2.5, 1.1, 1.0, 1.2, 3.2, 2.7)
B <- c(1.5, 2.7, 6.6, 4.6, 1.1, 1.2, 5.7, 3.2, 1.2, 1.3)
plot(density(A))    
lines(density(B))


## @knitr 
wilcox.test(A,B)


## @knitr 
xbar1 <- 79; xbar2 <- 110
n1 <- n2 <- 250
s1 <- 25; s2 <- 20
sp <- sqrt(( (n1-1)*s1^2 + (n2-1)*s2^2 )/(n1+n2-2))
SE <- sp * sqrt(1/n1 + 1/n2)
T <- (xbar1 - xbar2)/SE
2 * pt(T, df = n1+n2-2)       # T is negative, two sided


## @knitr 
xbar1 <- 5.3; xbar2 <- 5.4
sp <- 2.5                      # clearly as pooled value is average
n1 <- 200; n2 <- 207
SE <- sp*sqrt(1/n1 + 1/n2)
T <- (xbar1 - xbar2)/SE
2*pt(T, df = n1 + n2 - 2)     # T is negative


## @knitr fig.keep="none"
plot(age ~ dage, data=babies)
b = subset(babies, subset= dage < 99 & age < 99 )
plot(age ~ dage, data=b)


## @knitr 
with(b, t.test(dage - age, alt="greater"))


## @knitr 
plot(temperature ~ factor(gender), data=normtemp)
t.test(temperature ~ factor(gender), data=normtemp)


## @knitr 
pre <- c(17,12,20,12,20,21,23,10,15,17,18,18)
post <- c(19,25,18,18,26,19,27,14,20,22,16,18)
t.test(post - pre, mu = 0, alt="greater")


## @knitr 
method1 <- c(45.9,57.6,54.9,38.7,35.7,39.2,45.9,43.2,45.4,54.8)
method2 <- c(48.2,64.2,56.8,47.2,43.7,45.7,53.0,52.0,45.1,57.5)
t.test(method1 - method2, mu = 0, alt="two.sided")


## @knitr fig.keep="none"
library(MASS)                 # load data set 
names(shoes)
plot(B ~ A, data=shoes)       # related (not shown)
with(shoes, t.test(A,B,paired=TRUE)$p.value)
with(shoes, t.test(A,B)$p.value)


## @knitr 
with(Galton, t.test(child,parent, paired=TRUE))


## @knitr 
x <- c(284, 279, 289, 292, 287, 295, 285, 279, 306, 298)  
y <- c(298, 307, 297, 279, 291, 335, 299, 300, 306, 291)
var.test(x,y)


