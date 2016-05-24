
## @knitr echo=FALSE, out.width=singlewide
## Simulation of samples from a parent population, summary of sample means
set.seed(100)
gray <- "gray70"
g <- .35
red <- rgb(g, g, g, alpha=0.75)
max_y <- dnorm(0)

plot.new()
plot.window(xlim=c(-3,3), ylim=c(0, max_y))
axis(1)
title(main="10 samples of size 16")

x <- seq(-3, 3, length=500)
points(x, dnorm(x), main="XXX", ylab="", type="l", bty="L")

lines(c(0,0), c(0, max_y), lty=2) ## col=gray
lines(c(0,1), c(dnorm(1), dnorm(1)), lty=2) ##, col=gray)


n <- 16
m <- 10

ys <- max_y * seq(1/4, 3/4, length=m)

xbars <- sapply(ys, function(y) {
  x <- rnorm(n)
  points(x, y + 0*x, col=gray, pch=16)
  points(mean(x), y, col=red, pch=15, cex=2)
  mean(x)
})


## boxplot
## xmed <- mean(xbars)
## xiqr <- IQR(xbars)

## rect(xmed - xiqr/2, (1/8 - 1/16) * max_y, xmed + xiqr/2, (1/8 + 1/16) * max_y, border=red, cex=2)
## lines(xmed * c(1,1), c((1/8 - 1/16) * max_y, (1/8 + 1/16) * max_y), col=red, cex=2)
## lines(c(min(xbars), xmed - xiqr/2), 1/8 * max_y * c(1,1), col=red, cex=2)
## lines(c(max(xbars), xmed + xiqr/2), 1/8 * max_y * c(1,1), col=red, cex=2)


d <- density(xbars)
y <- d$y
lines(d$x, max_y *(1/8 - 2/16 + d$y/max(d$y)*3/16), col=red)

points(xbars, max_y*(1/8 - 2/16) + 0*xbars, cex=0.5, col=red)


## @knitr 
mu <- 100; sigma <- 16                  # population parameters
x <- rnorm(16, mean=mu, sd=sigma)       # our sample
mean(x)                                 # xbar


## @knitr echo=FALSE, out.width=singlewide
## create "wild" boxplot of simulation. Named after Chris Wild
d = stack(list("n=50"=rnorm(50), "n=250"=rnorm(250), "n=1000"=rnorm(1000)))
boxplot(values ~ ind, d, horizontal=TRUE, lwd=2)

QT = replicate(50, {
  d = stack(list("n=50"=rnorm(50), "n=250"=rnorm(250), "n=1000"=rnorm(1000)))
  boxplot(values ~ ind, d, horizontal=TRUE, border=rgb(.75, .75, .75, .75), add=TRUE)
})
boxplot(values ~ ind, d, horizontal=TRUE, lwd=2, add=TRUE)


## @knitr 
mu <- 100; sigma <- 16
M <- 4; n <- 16                         # change M for more
res <- numeric(M)
for (i in 1:M) {
  res[i] <- mean(rnorm(n, mean=mu, sd=sigma))
}
res


## @knitr 
xbar <- function(i)  mean(rnorm(n, mean=mu, sd=sigma)) #  function
sapply(1:M, xbar)


## @knitr 
Xbar <- Vectorize(xbar)
Xbar(1:M)


## @knitr 
replicate(M, mean(rnorm(n, mean=mu, sd=sigma))) # replicate(n, expr)


## @knitr 
x <- matrix(rnorm(M*n, mean=mu, sd=sigma), nrow=n) 
dim(x)                                  # M columns, n rows
apply(x, 2, mean)


## @knitr 
zstat <- function(x, mu, sigma) {
  (mean(x) - mu) / (sigma/sqrt(length(x)))
}


## @knitr zstat-normal-pop, eval=FALSE
## M <- 2000; n <- 7
## mu <- 100; sigma <- 16
## res <- replicate(M, {
##   x <- rnorm(n, mean=mu, sd=sigma)
##   zstat(x, mu, sigma)
## })
## qqnorm(res, main="Normal, n=7")


## @knitr zstat-exp-pop-n-7, eval=FALSE
## M <- 2000; n <- 7
## rate <- 2; mu <- sigma <- 1/rate
## res <- replicate(M, {
##   x <- rexp(n, rate=rate)
##   zstat(x, mu, sigma)
## })
## qqnorm(res, main="Exponential, n=7")


## @knitr zstat-exp-pop-n-150, eval=FALSE
## n <- 150
## res <- replicate(M, {
##   x <- rexp(n, rate=rate)
##   zstat(x, mu, sigma)
## })
## qqnorm(res, main="Exponential, n=150")


## @knitr echo=FALSE, out.width=triplewide
M <- 2000; n <- 7
mu <- 100; sigma <- 16
res <- replicate(M, {
  x <- rnorm(n, mean=mu, sd=sigma)
  zstat(x, mu, sigma)
})
qqnorm(res, main="Normal, n=7")
M <- 2000; n <- 7
rate <- 2; mu <- sigma <- 1/rate
res <- replicate(M, {
  x <- rexp(n, rate=rate)
  zstat(x, mu, sigma)
})
qqnorm(res, main="Exponential, n=7")
n <- 150
res <- replicate(M, {
  x <- rexp(n, rate=rate)
  zstat(x, mu, sigma)
})
qqnorm(res, main="Exponential, n=150")


## @knitr t-stat
tstat <- function(x, mu) (mean(x) - mu) / (sd(x) / sqrt(length(x)))


## @knitr eval=FALSE
## mu <- 0; sigma <- 1                     # use defaults
## M <- 750; n <- 4
## res <- replicate(M, tstat(rnorm(n, mu, sigma), mu))
## boxplot(res)


## @knitr echo=FALSE, out.width=singlewide
set.seed(1400)
ns <- c(4, 10, 20, 30, 50, 100)
res <- lapply(ns, function(n) replicate(M, tstat(rnorm(n, mu, sigma), mu)))
names(res) <- paste("n=", ns)
res[["normal"]] <- rnorm(M)
boxplot(res)


## @knitr mean-median-normal, eval=FALSE
## M <- 1000; n <- 35
## res_mean <- replicate(M, mean(rnorm(n)))
## res_median <- replicate(M, median(rnorm(n)))
## boxplot(list("sample mean"=res_mean, "sample median"=res_median),
##         main="Normal population")


## @knitr mean-median-exponential, eval=FALSE
## M <- 1000; n <- 35
## res_mean <- replicate(M, mean(rexp(n))) - 1
## res_median <- replicate(M, median(rexp(n))) - log(2)
## boxplot(list("sample mean"=res_mean, "sample median"=res_median),
##         main="Exponential population")


## @knitr echo=FALSE, out.width=doublewide
set.seed(10)
M <- 1000; n <- 35
res_mean <- replicate(M, mean(rnorm(n)))
res_median <- replicate(M, median(rnorm(n)))
boxplot(list("sample mean"=res_mean, "sample median"=res_median),
        main="Normal population")
M <- 1000; n <- 35
res_mean <- replicate(M, mean(rexp(n))) - 1
res_median <- replicate(M, median(rexp(n))) - log(2)
boxplot(list("sample mean"=res_mean, "sample median"=res_median),
        main="Exponential population")


## @knitr 
N <- 80; k <- 20; j <- 10
x <- sample(1:N, j, replace=FALSE)
x


## @knitr 
sum(x %in% 1:k)


## @knitr 
res <- replicate(10000, {
  x <- sample(1:N, j, replace=FALSE)
  sum(x %in% 1:k)
})


## @knitr 
c("zero matches"= sum(res == 0)/length(res), # prop.table to see all
  "one match"   = sum(res == 1)/length(res))


## @knitr 
sum(res >= 5) / length(res)             # 5 or more matches


## @knitr 
res <- replicate(1000, sum(sample(1:6, 2, replace=TRUE)))
sum(res == 7) > sum(res == 8)


## @knitr 
sam <- sample(rivers, 10, replace=TRUE)


## @knitr 
sams <- replicate(1000,  mean(sample(rivers, 10, replace=TRUE)))


## @knitr 
c("mean of samples"=mean(sams), mean=mean(rivers))


## @knitr 
ctrl <- c(23, 33, 40)
treatment <- c(19, 22, 25, 26)
the_data <- stack(list(ctrl=ctrl, treatment=treatment))
aggregate(values ~ ind, the_data, mean)


## @knitr 
cmbs <- combn(7, 3)                     # 35 possibilities
cmbs[, 1:6]                             # first 6 columns


## @knitr 
i <- 1
ind <- cmbs[,1]                         # in columns
obs <- mean(the_data$value[ind]) - mean(the_data$value[-ind])
obs


## @knitr 
res <- apply(cmbs, 2, function(ind) {
  mean(the_data$value[ind]) - mean(the_data$value[-ind])
})


## @knitr 
sum(res >= obs)
sum(res >= obs) / length(res)


## @knitr 
caf <- c(245, 246, 246, 248, 248, 248, 250, 250, 250, 252)
no_caf <- c(242, 242, 242, 244, 244, 245, 246, 247, 248, 248)
the_data <- stack(list(caffeine=caf, no_caffeine=no_caf))


## @knitr 
obs <- mean(caf) - mean(no_caf)
obs


## @knitr 
sample(1:20, 10, replace=FALSE)


## @knitr echo=FALSE
set.seed(314)


## @knitr 
res <- replicate(2000, {
  ind <- sample(1:20, 10, replace=FALSE)
  mean(the_data$value[ind]) - mean(the_data$value[-ind])
})


## @knitr 
sum(res > obs) / length(res)


## @knitr 
mu <- 100; sigma <- 16
M <- 1000; n <- 4
res <- replicate(4, mean(rnorm(n, mu, sigma)) - mu)


## @knitr echo=FALSE
set.seed(10)


## @knitr 
mu <- 100; sigma <- 16
M <- 1000; n <- 4

res <- replicate(M, {
  x <- rnorm(n, mu, sigma)
  SE <- sd(x)/sqrt(n)                   # standard error
  (mean(x) - mu) / SE
})


## @knitr 
quantile(res, c(0.025, 0.975))


## @knitr 
diabetes <- subset(Medicare, 
               subset= DRG.Definition =="638 - DIABETES W CC")
gap <- with(diabetes, Average.Covered.Charges-Average.Total.Payments)


## @knitr 
range(gap)


## @knitr 
xbar <- mean(gap)
xbar


## @knitr 
xstar <- sample(gap, length(gap), replace=TRUE)
head(xstar)


## @knitr 
mean(xstar) - xbar


## @knitr 
M <- 2000
res <- replicate(M, {
  xstar <- sample(gap, length(gap), replace=TRUE)
  mean(xstar) - xbar
})


## @knitr 
alpha <- 0.05
xbar + quantile(res, c(alpha/2, 1-alpha/2))


## @knitr prior
p <- seq(0.05, 0.95, by=0.1)            # possible values
prior <- c(2, 4, 8, 8, 4, 2, 1, 1, 1, 1) # how likely
prior <- prior / sum(prior)             # make a probability
prior


## @knitr 
like <- p^11 * (1-p)^(27 - 11)


## @knitr 
posterior <- like * prior
posterior <- posterior / sum(posterior)


## @knitr prior-posterior, echo=FALSE, out.width=doublewide
spike_plot <- function(p, x, ...) {
  plot(p, x, type="h", ...)
  points(p, x, pch=16, cex=1.5)
}
spike_plot(p, prior, ylab="Prior distribution")
spike_plot(p, posterior, ylab="Posterior distribution")


## @knitr 
ind <- 0 <= p & p <= 0.5
sum(prior[ind])


## @knitr 
sum(posterior[ind])


