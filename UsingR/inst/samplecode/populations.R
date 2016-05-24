##' plot function f over a,b with shading in [l,r]
shade <- function(f, a, b, l=a, r=b, plot=TRUE, col="gray60", ...) {
  if(plot)
    curve(f, a, b, ...)

  x <- seq(l, r, length=250)
  y <- f(x)
  polygon(c(x, rev(x)), c(y, 0*x), col=col, ...)
}




## @knitr 
tbl <- xtabs(~ Sex + Smoke, data=survey)
tbl


## @knitr 
margin.table(tbl, margin=1)
118 / (118 + 117)                       # sum(tbl[1,]) / sum(tbl)


## @knitr spike-plot, eval=FALSE
## k <- 0:4
## p <- c(1, 2, 3, 2, 1); p <- p/sum(p)
## plot(k, p, type="h", xlab="k", ylab="probability", ylim=c(0,max(p)))
## points(k, p, pch=16, cex=2)      # add the balls to top of spike


## @knitr echo=FALSE, out.width=doublewide
k <- 0:4
p <- c(1, 2, 3, 2, 1); p <- p/sum(p)
plot(k, p, type="h", xlab="k", ylab="probability", ylim=c(0,max(p)))
points(k, p, pch=16, cex=2)      # add the balls to top of spike
## spinner
op <- par(no.readonly = TRUE)
par(mai=c(0,0,0,0))
pie(p, col="gray80")
r <- 0.5
t <- 2*pi*runif(1)
arrows(r * cos(t), r*sin(t), -r*cos(t), -r*sin(t), lwd=3)
par(op)


## @knitr echo=FALSE
set.seed(4)


## @knitr 
k <- 0:2
p <- c(1, 2, 1); p <- p/sum(p)          # add to 1
sample(k, size=1, prob=p)               
sample(k, size=1, prob=p)               # likely different


## @knitr 
sample(1:6, size=1) + sample(1:6, size=1) 


## @knitr echo=FALSE, out.width=singlewide
## drawa simple shade graphic
plot.new()
plot.window(xlim=c(0,15), ylim=c(0,.25))
f <- function(x) dchisq(x, df=3)
curve(f, add=TRUE)
b <- 10
x <- seq(0, b, length=250)
polygon(c(x, rev(x)), c(f(x), 0*x), col="gray60")
abline(h=0)

text(b, 0, "b", pos=1)
text(2, 0.075, expression(P(X <= b)))


## @knitr echo=FALSE, out.height=".33\\linewidth", out.width="\\linewidth"
plot.new()
plot.window(xlim=c(0, 3*15-1), ylim=c(-.05,.25))
f <- function(x) dchisq(x, df=3)
col = "gray60"

## plots f's
x <- seq(0, 14, length=250)
lines(x, f(x))
lines(x + 15, f(x))
lines(x + 30, f(x))

lines(c(0,14), c(0,0))
lines(15+c(0,14), c(0,0))
lines(30+c(0,14), c(0,0))

## shade
a = 2; b=10

x = seq(a, b, length=250)
polygon(c(x, rev(x)), c(f(x), 0*x), col=col)

x = seq(0, b, length=250)
polygon(15 + c(x, rev(x)), c(f(x), 0*x), col=col)

x = seq(0, a, length=250)
polygon(30 + c(x, rev(x)), c(f(x), 0*x), col=col)

## label
text(c(a,b), c(0,0),      c("a", "b"), pos=1)
text(15 + c(a,b), c(0,0), c("a", "b"), pos=1)
text(30 + c(a,b), c(0,0), c("a", "b"), pos=1)

## math
text(15 -2, 0.25/2, expression("="), cex=1.5, lwd=2)
text(30 -2, 0.25/2, expression("-"), cex=1.5, lwd=2)



## @knitr 
## toss a coin 10 times. Heads=1, tails=0
sample(0:1,size=10,replace=TRUE) 
sample(1:6,size=10,replace=TRUE) ## roll a die 10 times
## sum of roll of a pair of dice roll 10 times
sample(1:6,size=10,replace=TRUE) + sample(1:6,size=10,replace=TRUE)


## @knitr 
sample(rep(0:1, c(3200,6800)), size=10, replace=TRUE)


## @knitr 
p <- 0.68
sample(0:1, size=10, replace=TRUE, prob=c(1-p, p))


## @knitr 
 max(sample(1:6,1),sample(1:6,1))


## @knitr 
k <- 1:6
sample(k, 1, prob= (2*k-1)/36)


## @knitr 
with(nba.draft, sample(Team, 1, prob=Balls))


## @knitr 
dunif(x=1, min=0, max=3)            # 1/3 of area is to left of 1
punif(q=2, min=0, max=3)            # 1/(b-a) is 2/3
qunif(p=1/2, min=0, max=3)          # half way between 0 and 3
runif(n=1, min=0, max=3)            # a random value in [0,3]


## @knitr 
ps <- seq(0, 1, by=.2)                        # probabilities
names(ps) <- as.character(seq(0, 100, by=20)) # give names
qunif(ps, min=0, max=1)                       


## @knitr 
runif(5, min=0, max=1:5)                # recycles min


## @knitr relate-d-and-r, eval=FALSE
## x <- runif(100)                         # large sample, 1000 points
## d <- density(x)
## curve(dunif, -0.1, 1.1,
##       ylim=c(0, max(d$y, 1)))           # plots function
## lines(d, lty=2)                         # add density estimate
## rug(x)                                  # indicates sample


## @knitr echo=FALSE, out.width=singlewide
set.seed(18)
x <- runif(100)                         # large sample, 1000 points
d <- density(x)
curve(dunif, -0.1, 1.1, 
      ylim=c(0, max(d$y, 1)))           # plots function
lines(d, lty=2)                         # add density estimate
rug(x)                                  # indicates sample


## @knitr 
n <- 10; p <- 1/4
sample(0:1, size=n, replace=TRUE, prob=c(1-p, p))


## @knitr 
choose(10,5) * (1/2)^5 * (1/2)^(10-5)


## @knitr 
dbinom(5, size=10, prob=1/2)    


## @knitr 
sum(dbinom(0:6, size=10, prob=1/2))
pbinom(6, size=10, p=1/2)


## @knitr 
sum(dbinom(7:10,size=10,prob=1/2))
1 - pbinom(6,size=10,p=1/2)
pbinom(6,size=10,p=1/2, lower.tail=FALSE) # k = 6 not 7!


## @knitr binomial-spike-plot, eval=FALSE
## n <- 10; p <- 1/2
## heights <- dbinom(0:10, size=n, prob=p)
## plot(0:10, heights, type="h",
##  main="Spike plot of X", xlab="k", ylab="p.d.f.")
## points(0:10, heights, pch=16, cex=2)


## @knitr echo=FALSE, out.width=doublewide
n <- 10; p <- 1/2
heights <- dbinom(0:10, size=n, prob=p)
plot(0:10, heights, type="h",
 main="Spike plot of X", xlab="k", ylab="p.d.f.")    
points(0:10, heights, pch=16, cex=2)
  n <- 10; p <- 1/2
ps <- pbinom(0:n, size=n, p=p)
plot(stepfun(0:10, c(0, ps)), 
     verticals=FALSE,
     main="c.d.f. of binomial(n=10, p=1/2)"
     )


## @knitr 
pbinom(60, size=100, prob=0.62)


## @knitr 
pnorm(1, mean=0, sd=1)      
pnorm(4.5, mean=4, sd=1/2)              # same z-score as above


## @knitr echo=FALSE, out.width=doublewide
f <- function(mu, sd) function(x) dnorm(x, mean=mu, sd=sd)
a <- -3; b <- 5.5
mu <- 0; sd <- 1
x <- seq(a, b, length=500)
plot(x, f(mu, sd)(x), bty="L", type="l", 
     ylim=c(0, dnorm(0, sd=1/2)),
     ylab="", xlab="")

x <- seq(a, mu + 3/2*sd, length=250)
polygon(c(x, rev(x)), c(f(mu,sd)(x), 0*x), col="gray60")

lines(c(mu,mu), c(0, dnorm(mu, mean=mu, sd=sd)), lty=2)
lines(c(mu, mu + sd), dnorm(mu+sd, mean=mu, sd=sd)*c(1,1), lty=2)
text(mu + 1/2*sd, dnorm(mu+sd, mean=mu, sd=sd), as.character(sd), pos=1)


mu <- 4; sd <- 1/2
x <- seq(a, b, length=500)
lines(x, f(mu, sd)(x))
x <- seq(a, mu + 3/2*sd, length=250)
polygon(c(x, rev(x)), c(f(mu, sd)(x), 0*x), col="gray60")


lines(c(mu,mu), c(0, dnorm(mu, mean=mu, sd=sd)), lty=2)
lines(c(mu, mu + sd), dnorm(mu+sd, mean=mu, sd=sd)*c(1,1), lty=2)
text(mu + 1/2*sd, dnorm(mu+sd, mean=mu, sd=sd), as.character(sd), pos=1)




## @knitr 
qnorm(c(0.25, 0.5, 0.75))


## @knitr 
pnorm(1) - pnorm(-1)


## @knitr 
1 - 2*pnorm(-2)               # subtract area of two tails
diff(pnorm(c(-3, 3)))          # use diff to subtract


## @knitr 
mu <- 70.2; sigma <- 2.89   
pnorm(72, mean=mu, sd=sigma)


## @knitr 
conv <- 0.0254
pnorm(2/conv, mean=mu, sd=sigma)


## @knitr 
p = 1 - 1 / (3.5e9)
qnorm(p, mu, sigma)/12


## @knitr echo=FALSE
set.seed(10)


## @knitr 
mu <- 100; sigma <- 10
res <- rnorm(1000, mean=mu, sd=sigma)


## @knitr 
k <- 1; sum(res > mu - k*sigma & res < mu + k*sigma)
k <- 2; sum(res > mu - k*sigma & res < mu + k*sigma)
k <- 3; sum(res > mu - k*sigma & res < mu + k*sigma)


## @knitr hist-uniform-boxplot, eval=FALSE
## res = runif(50, min=0, max=10)
## ## fig= setting uses bottom 35% of diagram
## par(fig=c(0,1,0,.35))
## boxplot(res,horizontal=TRUE, bty="n", xlab="uniform sample")
## ## fig= setting uses top 75% of figure
## par(fig=c(0,1,.25,1), new=TRUE)
## hist(res, prob=TRUE, main="", col=gray(.9))
## lines(density(res),lty=2)
## curve(dunif(x, min=0, max=10), lwd=2, add=TRUE)
## rug(res)


## @knitr hist-exp-boxplot, echo=FALSE, eval=FALSE
## res = rexp(50, 5)
## ## fig= setting uses bottom 35% of diagram
## par(fig=c(0,1,0,.35))
## boxplot(res,horizontal=TRUE, bty="n", xlab="exponential sample")
## ## fig= setting uses top 75% of figure
## par(fig=c(0,1,.25,1), new=TRUE)
## hist(res, prob=TRUE, main="", col=gray(.9))
## lines(density(res),lty=2)
## curve(dexp(x, 5), lwd=2, add=TRUE)
## rug(res)


## @knitr echo=FALSE, out.width=doublewide
res = runif(50, min=0, max=10)
## fig= setting uses bottom 35% of diagram
par(fig=c(0,1,0,.35))
boxplot(res,horizontal=TRUE, bty="n", xlab="uniform sample")
## fig= setting uses top 75% of figure
par(fig=c(0,1,.25,1), new=TRUE)
hist(res, prob=TRUE, main="", col=gray(.9))
lines(density(res),lty=2)
curve(dunif(x, min=0, max=10), lwd=2, add=TRUE)
rug(res)
res = rexp(50, 5)
## fig= setting uses bottom 35% of diagram
par(fig=c(0,1,0,.35))
boxplot(res,horizontal=TRUE, bty="n", xlab="exponential sample")
## fig= setting uses top 75% of figure
par(fig=c(0,1,.25,1), new=TRUE)
hist(res, prob=TRUE, main="", col=gray(.9))
lines(density(res),lty=2)
curve(dexp(x, 5), lwd=2, add=TRUE)
rug(res)


## @knitr echo=FALSE, out.width=singlewide
set.seed(1)
res = rlnorm(50, 0, 1)
## fig= setting uses bottom 35% of diagram
par(fig=c(0,1,0,.35))
boxplot(res,horizontal=TRUE, bty="n", xlab="Log-normal sample")
## fig= setting uses top 75% of figure
par(fig=c(0,1,.25,1), new=TRUE)
hist(res, prob=TRUE, main="", col=gray(.9))
lines(density(res),lty=2)
curve(dlnorm(x, 0, 1), lwd=2, add=TRUE)
rug(res)


## @knitr 
qt(c(0.025, 0.975), df=10)              # 10 degrees of freedom
qf(c(0.025, 0.975), df1=10, df2=5)      # 10 and 5 degrees of freedom
qchisq(c(0.025, 0.975), df=10)          # 10 degrees of freedom


## @knitr 
sum(dbinom(3:5, size=5, prob=1/2))


## @knitr 
dbinom(12, size=12, prob=.3)


## @knitr 
sum(dbinom(k, size=1e6, prob=1/2))


## @knitr 
dbinom(4, size=4, prob=1/3)


## @knitr 
1 - pbinom(0,4,1/6)           # one or more 6 in 4 rolls
1 - pbinom(0,24,1/36)         # one or more double sixes


## @knitr 
n <- 18:30                              # some range, may be wrong
names(n) <- n                           # helps, not necessary
pbinom(0, size=n, p=1/36)               # use vectorized n


## @knitr 
pbinom(35, size=100, prob=0.4)      


## @knitr 
pnorm(2.2)
pnorm(2) - pnorm(-1)
pnorm(2.5, lower.tail=FALSE)  # or 1 - pnorm(2.5)
qnorm(.95)


## @knitr 
1 - pnorm(450, mean=350, sd=75)


## @knitr 
1 - pnorm(22, mean=20.6, sd=5.5)


## @knitr 
1000000 * diff(pnorm(c(22,23), mean=20.6, sd=5.5))      


## @knitr 
pnorm(26, mean=24.9, sd=1.05)
qnorm(0.85, mean=24.9, sd=1.05)


## @knitr 
mu <- 3.20
sigma <- 0.35
pnorm(4, mean=mu, sd=sigma) - pnorm(3.5, mean=mu, sd=sigma)


## @knitr 
pnorm(6) - pnorm(-6)
1 - (pnorm(6) - pnorm(-6))


## @knitr 
x <- father.son$fheight
x <- (x - mean(x))/sd(x)                # or scale(x)[,1]
sum(abs(x) < 1)/length(x)
sum(abs(x) < 2)/length(x)
sum(abs(x) < 3)/length(x)


## @knitr 
ps <- seq(0, 1, by=0.2)
qnorm(ps)


## @knitr 
mu <- 1/2
sigma <- sqrt(1/12)
k <- 1; diff(punif(mu + c(-1,1)*k*sigma))
k <- 2; diff(punif(mu + c(-1,1)*k*sigma))
k <- 3; diff(punif(mu + c(-1,1)*k*sigma))


## @knitr 
mu <- 5
sigma <- 5
k <- 1; diff(pexp(mu + c(-1,1)*k*sigma, rate=1/mu))
k <- 2; diff(pexp(mu + c(-1,1)*k*sigma, rate=1/mu))
k <- 3; diff(pexp(mu + c(-1,1)*k*sigma, rate=1/mu))


## @knitr eval=FALSE
## qqnorm(runif(100))
## qqnorm(rt(100, df=3))
## qqnorm(rnorm(100))


## @knitr eval=FALSE
## curve(dnorm(x), -4, 4)


## @knitr eval=FALSE
## k <- 5; curve(dt(x, df=k), lty=k, add=TRUE)


## @knitr eval=FALSE
## curve(dnorm(x),-4,4)
## k <-  5;   curve(dt(x,df=k), lty=k, add=TRUE)
## k <- 10;  curve(dt(x,df=k), lty=k, add=TRUE)
## k <- 25;  curve(dt(x,df=k), lty=k, add=TRUE)
## k <- 50;  curve(dt(x,df=k), lty=k, add=TRUE)
## k <- 100; curve(dt(x,df=k), lty=k, add=TRUE)


## @knitr eval=FALSE
## curve(dchisq(x,df=2), 0, 100)


## @knitr eval=FALSE
## k <- 8; curve(dchisq(x,df=k), add=TRUE)


## @knitr normal-clt, eval=FALSE
## n <- 25; curve(dnorm(x, mean=0, sd=1/sqrt(n)), -3, 3,
##                xlab="x", ylab="Densities of sample mean", bty="l")
## n <- 5;  curve(dnorm(x, mean=0, sd=1/sqrt(n)), add=TRUE)
## n <- 1;  curve(dnorm(x, mean=0, sd=1/sqrt(n)), add=TRUE)


## @knitr echo=FALSE, out.width=singlewide
n <- 25; curve(dnorm(x, mean=0, sd=1/sqrt(n)), -3, 3, 
               xlab="x", ylab="Densities of sample mean", bty="l")
n <- 5;  curve(dnorm(x, mean=0, sd=1/sqrt(n)), add=TRUE)
n <- 1;  curve(dnorm(x, mean=0, sd=1/sqrt(n)), add=TRUE)


## @knitr 
mu <- 70.2; sigma <- 2.89; n <- 25
diff(pnorm(70:71, mu, sigma/sqrt(n)))


## @knitr 
diff( pnorm(70:71, mu, sigma) )


## @knitr echo=FALSE, out.width=singlewide
n <- c( 5,25,100)
x <- sapply(n, function(i) replicate(2000,mean(rexp(i))))
l <- apply(x, 2, density)
xlim <- range(sapply(l, function(d) range(d$x)))
ylim <- range(sapply(l, function(d) range(d$y)))
plot(l[[1]], xlim=xlim, ylim=ylim, main="Density estimates of exponentials")
QT <- sapply(l[-1], lines)


## @knitr 
pnorm(0.9, mean=1, sd = 1/sqrt(20))


## @knitr echo=FALSE, out.width=singlewide
## plot binomial and normal curve to show the approximation

k <- 10:28;
n <- 30; 
p <- 2/3

mu <- n*p
sigma <- sqrt(n*p*(1-p))

heights <- dbinom(k,n,p)

plot(k,heights,type="h",
     bty="L",
     xlab = "binomial(30,2/3)",
     ylab = "probability"
     )

## add boxes
plot.box = function(x,y,col=jvgray) {
  x = c(x-1/2,x-1/2,x+1/2,x+1/2)
  y = c(0,y,y,0)
  polygon(x,y,
          lty=2,
#          density = 4,
          col=col
          )
}

for(i in 1:length(k)) {
  col = gray(.8)
  if(k[i] >= 23) {
    col = gray(1)
  }
  plot.box(k[i],heights[i],col=col)
}



points(k,heights,pch=16,cex=2)

## add curve
f <- function(x) {
  dnorm(x,mu,sigma)
}

curve(f,14,28,add=T)


## @knitr 
n <- 100; p <- 1/2; b <- 42
pbinom(b,n,p)
pnorm(b+1/2, n*p, sqrt(n*p*(1-p)))


## @knitr 
n <- 600; p <- .300
phat <- .350; a <-  n*phat
1 - pnorm(a - 1/2, n*p, sqrt(n*p*(1-p)))


## @knitr 
n <- 1000; p <- 1/2
mu <- n*p; sigma <- sqrt(n*p*(1-p))
1 - pnorm(550 - 1/2, mu, sigma)


## @knitr 
1 - pnorm(3500/15, mean=180, sd=25/sqrt(15))      


## @knitr 
1 - pnorm(775/30, mean=25, sd=2/sqrt(30))


## @knitr 
pnorm(75/21, mean=4, sd = 1/sqrt(21))


