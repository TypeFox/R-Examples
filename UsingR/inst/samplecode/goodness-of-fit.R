## @knitr 
cols <- c("blue", "brown", "green", "orange", "red", "yellow", 
          "purple")
prob <- c(1, 1, 1, 1, 2, 2, 2)       # ratio of colors
prob <- prob / sum(prob)
n <- 30
bagful <- sample(cols, n, replace=TRUE, prob=prob)
table(bagful)


## @knitr 
choose(100,35)*choose(65,40)*choose(25,25) * .35^35 * .35^40 * .30^25


## @knitr chi-sq-simulation, echo=FALSE, fig.out=singlewide
set.seed(100)
M <- 2000; n <- 20
p <- c(3,4,5)/12

res <- replicate(M, {
  x <- sample(1:3, n, replace=TRUE, prob=p)
  y <- sapply(1:3, function(i) sum(x==i))
  expected <- n * p
  chi <- sum( (y - expected)^2/expected )
  chi
})

col <- rgb(.7,.7,.7,.75)
hist(res, prob=TRUE, breaks=50, col=col, ylab="Probability", xlab="res", main="Chi-squared simulation")
curve(dchisq(x, df=length(p)-1), add=TRUE, lwd=2)
              



## @knitr 
y <- c(35, 40, 25)
p <- c(35, 35, 30)               # ratios
p <- p/sum(p)                    # proportions
n <- sum(y)
chi2 <- sum( (y - n*p)^2 / (n*p) )
chi2
pchisq(chi2, df=3 - 1, lower.tail=FALSE) 


## @knitr 
chisq.test(y, p=p)


## @knitr 
library(UsingR)  
amt <- with(samhda, amt.smoke[amt.smoke < 98])
y <- table(amt)
y
ps <- c(0.15, 0.05, 0.05, 0.05, 0.10, 0.20, 0.40)    
chisq.test(y, p=ps)


## @knitr 
y <- c(35, 40, 25)
n <- sum(y)
phat1 <- phat2 <- sum(y[1:2])/(2*n)
phat3 <- 1 - phat1 - phat2
phat <- c(phat1, phat2, phat3)
#
obs <- sum((y - n*phat)^2/(n*phat))
obs
pchisq(obs, df =1 , lower.tail=FALSE)


## @knitr 
x <- c(13, 17, 9, 17, 18, 26)
chisq.test(x)


## @knitr 
obs <- c(315, 197, 141, 39, 16, 79)
p <- c(.486, .315, .125, .028, .006, .040)
chisq.test(obs, p=p)


## @knitr 
prob <- mandms["milk chocolate", ]
prob <- unlist(prob)
bagful <- c(15, 34, 7, 19, 29, 24) 
names(bagful) = c("blue", "brown", "green", "orange", "red", "yellow")
chisq.test(bagful, p=prob/sum(prob))


## @knitr 
prob <- mandms["Peanut", ]
prob <- unlist(prob)
chisq.test(bagful, p=prob/sum(prob))


## @knitr 
chisq.test(table(pi2000))


## @knitr 
counts <- c(28,39,23,22,11)
freq <- c(9,12,9,8,4)
chisq.test(counts, p=freq/sum(freq))


## @knitr 
all.names <- paste(bright.stars$name, sep="", collapse="")
x <- unlist(strsplit(tolower(all.names), ""))
letter.dist <- sapply(letters, function(i) sum(x == i))


## @knitr 
ps <- scrabble$frequency[1:26]
ps <- ps/sum(ps) 


## @knitr 
all.names <- paste(bright.stars$name, sep="", collapse="")
x <- unlist(strsplit(tolower(all.names), ""))
letter.dist <- sapply(letters, function(i) sum(x == i))
ps <- scrabble$frequency[1:26]
chisq.test(letter.dist, p=ps/sum(ps))


## @knitr 
murder <- c(63 , 53 , 50 , 51 , 55 , 52 , 56)
chisq.test(murder)


## @knitr 
n <- sum(murder)
phatw <- (53 + 65)/(2*n)
phatd <- (1 - 2*phatw)/5
e <- n * c(phatw, rep(phatd,5), phatw)
cs <- sum ( (murder - e)^2/e )
cs
1 - pchisq(cs, df=1)


## @knitr 
colors <- c(41,48,105,58)
n <- sum(colors)
phat <- mean(41/n,48/n)
yellowhat <- 105/n;
greenhat <- 58/n
exp <- n*c(phat, phat, yellowhat, greenhat)
cs <- sum( (colors - exp)^2/exp )
cs
pchisq(cs,df=2,lower.tail=FALSE)


## @knitr fig.keep="none"
M <- 2000; n <- 20
p <- c(3,4,5)/12

res <- replicate(M, {
  x <- sample(1:3, n, replace=TRUE, prob=p)
  y <- sapply(1:3, function(i) sum(x==i))
  expected <- n * p
  chi <- sum( (y - expected)^2/expected )
  chi
})

col <- rgb(.7,.7,.7,.75)
hist(res, prob=TRUE, breaks=50, 
     col=col, ylab="Probability", xlab="res",
     main="Chi-squared simulation")
curve(dchisq(x, df=length(p)-1), add=TRUE, lwd=2)


## @knitr 
seatbelt <- rbind(c(56,8), c(2,16))
seatbelt
chisq.test(seatbelt)


## @knitr 
tbl <- xtabs( ~ gender + amt.smoke,    # no left side in formula
             subset = amt.smoke < 98 & gender !=7, 
             data=samhda)
tbl
chisq.test(tbl)


## @knitr 
chisq.test(tbl,simulate.p.value=TRUE)


## @knitr 
celexa <- c(2, 3, 7)
placebo <- c(2, 8, 2)
x <- rbind(celexa, placebo)
colnames(x) <- c("worse", "same", "better")
x
chisq.test(x)


## @knitr 
chisq.test(x, simulate.p.value=TRUE)


## @knitr 
accidents <- cbind(
  none=c(67,42,75,56,57),
  minor=c(10,6,8,4,15),
  major=c(5,5,4,6,1))
rownames(accidents) <- c("<18", "18-25", "26-40", "40-65", "65>")
accidents
chisq.test(accidents)


## @knitr 
aq <- airquality[complete.cases(airquality),]
aq <- transform(aq, 
  te = cut(Temp, quantile(Temp)),
  oz = cut(Ozone,quantile(Ozone))
)
xtabs(~ te + oz, data=aq)


## @knitr 
aq <- airquality[complete.cases(airquality),]
aq <- transform(aq, 
  te = cut(Temp, quantile(Temp)),
  oz = cut(Ozone,quantile(Ozone))
)
xtabs(~ te + oz, data=aq)
chisq.test(xtabs(~ te + oz, data=aq))


## @knitr 
sb.yes <- c(12813, 647, 359, 42)      
sb.no <- c(65963, 4000, 2642, 303)
chisq.test(rbind(sb.yes,sb.no))


## @knitr 
oral.lesion
chisq.test(oral.lesion)
chisq.test(oral.lesion, simulate.p.value=TRUE)


## @knitr 
retention <- rbind(
  nonblock=c(18, 15, 5, 8, 4),
  block = c(10, 5, 7, 18, 10))
colnames(retention) <- c(1:4, "5+")
retention
chisq.test(retention)


## @knitr 
y2011 <- c(63, 53, 50, 51, 55, 52, 56)
y2003 <- c(53, 42, 51, 45, 36, 37, 65)
x <- rbind(y2011, y2003)
chisq.test(x)


## @knitr echo=FALSE, out.width=doublewide
y <- rnorm(20)
plot(density(y), main="Densities")
curve(dnorm(x), add=TRUE, lty=2)
plot(ecdf(y), main="C.d.f.s")
curve(pnorm(x), add=TRUE, lty=2)


## @knitr echo=FALSE, out.width=singlewide
l <- list()
M <- 2000; n <- 25

set.seed(100)

l$normal <-  replicate(M, ks.test(rnorm(n), "pnorm")$statistic)
l$exponential <-  replicate(M, ks.test(rexp(n), "pexp")$statistic)
PT <- function(x) pt(x, df=5)
l$t <-  replicate(M, ks.test(rt(n, df=5), "PT")$statistic)
l$log_normal <- replicate(M, ks.test(rlnorm(n), "plnorm")$statistic)

d <- lapply(l, density)

plot(d[[1]], main="Distribution of KS statistic simulations", xlab="")
for (i in 2:length(d)) {
  lines(d[[i]], lty=i)
}

legend(.3, 8, legend=names(l), lty=1:5)



## @knitr 
x <- rnorm(100, mean=5, sd=2)
ks.test(x,"pnorm", mean=0, sd=2)          # "wrong" parameters
ks.test(x,"pnorm", mean=5, sd=2)$p.value  # correct parameters
x = runif(100, min=0, max=5)
ks.test(x,"punif", min=0, max=6)$p.value  # "wrong" parameters
ks.test(x,"punif", min=0, max=5)$p.value  # correct parameters


## @knitr 
library(UsingR)
sat.m <- stud.recs$sat.m; sat.v <- stud.recs$sat.v


## @knitr sat-eda, eval=FALSE
## boxplot(list(math=sat.m, verbal=sat.v), main="SAT scores")
## qqplot(sat.m, sat.v, main="Math and verbal SAT scores")
## plot(ecdf(sat.m), main="Math and verbal SAT scores")
## lines(ecdf(sat.v), lty=2)


## @knitr 
ks.test(sat.m,sat.v)


## @knitr echo=FALSE, out.width=triplewide
boxplot(list(math=sat.m, verbal=sat.v), main="SAT scores")
qqplot(sat.m, sat.v, main="Math and verbal SAT scores")
plot(ecdf(sat.m), main="Math and verbal SAT scores")
lines(ecdf(sat.v), lty=2)


## @knitr 
set.seed(100)


## @knitr ks-estimated-parameters, eval=FALSE
## res <- replicate(2000, {
##   x <- rnorm(25, mean=0, sd=1)
##   c(ks.test(x,pnorm, mean=mean(x), sd=sd(x))$statistic,
##     ks.test(x,pnorm, mean=0,       sd=1)$statistic)
## })
## plot(density(res[1,]), main="K-S sampling distribution", ylab="")
## lines(density(res[2,]), lty=2)
## legend(0.2, 12, legend=c("estimated", "exact"), lty=1:2)
## 


## @knitr echo=FALSE, out.width=singlewide
res <- replicate(2000, {
  x <- rnorm(25, mean=0, sd=1)
  c(ks.test(x,pnorm, mean=mean(x), sd=sd(x))$statistic,
    ks.test(x,pnorm, mean=0,       sd=1)$statistic)
})
plot(density(res[1,]), main="K-S sampling distribution", ylab="")
lines(density(res[2,]), lty=2)
legend(0.2, 12, legend=c("estimated", "exact"), lty=1:2)
                  


## @knitr 
shapiro.test(stud.recs$sat.m)
shapiro.test(stud.recs$sat.v)


## @knitr 
shapiro.test(OBP)$p.value


## @knitr 
shapiro.test(OBP[OBP<.5])$p.value


## @knitr 
library(MASS)
fitdistr(babyboom$wt,"normal")


## @knitr 
inter = diff(babyboom$running.time)


## @knitr 
out <- fitdistr(inter,"gamma")
out


## @knitr babyboom-inter-arrival-times, eval=FALSE
## plot(density(inter), ylim=c(0,0.025),
##      main="Compare estimated densities", xlab="inter")
## curve(dgamma(x, shape=out$estimate["shape"],
##              rate=out$estimate["rate"]), add=TRUE, lty=2)
## legend(100,.020,legend=c("density()","fitdistr()"),lty=1:2)
## #
## plot(ecdf(inter),
##  main="Compare ecdf with estimated cdf", xlab="inter")
## curve(pgamma(x,shape=1.208593, rate=0.036350), add=TRUE)
## legend(70,.8,legend=c("ecdf","estimated cdf"),lty=1:2)


## @knitr echo=FALSE, out.width=doublewide
plot(density(inter), ylim=c(0,0.025), 
     main="Compare estimated densities", xlab="inter")
curve(dgamma(x, shape=out$estimate["shape"], 
             rate=out$estimate["rate"]), add=TRUE, lty=2)
legend(100,.020,legend=c("density()","fitdistr()"),lty=1:2)
#
plot(ecdf(inter), 
 main="Compare ecdf with estimated cdf", xlab="inter")
curve(pgamma(x,shape=1.208593, rate=0.036350), add=TRUE)
legend(70,.8,legend=c("ecdf","estimated cdf"),lty=1:2)


## @knitr 
shapiro.test(babies$ht[babies$ht < 99])$p.value
shapiro.test(babies$wt[babies$wt < 999])$p.value


## @knitr fig.keep="none"
hist(brightness, prob=TRUE)
lines(density(brightness))
curve(dnorm(x, mean(brightness), sd(brightness)), add=TRUE)    


## @knitr 
shapiro.test(brightness)


## @knitr 
shapiro.test(normtemp$temperature)


## @knitr fig.keep="none"
library(MASS)
fitdistr(rivers,"gamma")
plot(density(rivers))
curve(dgamma(x,shape=2.578, rate= 0.00436), add=TRUE, lty=2)


## @knitr fig.keep="none"
qqplot(rivers, rgamma(100, shape=2.578, rate= 0.00436))      


## @knitr 
fitdistr(stud.recs$sat.m, "normal")


## @knitr 
res <- replicate(1000,  ks.test(rt(25,df=3),"pnorm")$p.value)


## @knitr echo=FALSE
set.seed(100)


## @knitr 
M <- 1000
res.t   <- replicate(M, ks.test(rt(25,df=3), "pnorm")$p.value)
res.exp <- replicate(M, ks.test(rexp(25)-1,  "pnorm")$p.value)
sum(res.t < .05)/M
sum(res.exp < .05)/M


## @knitr fig.keep="none"
qqplot(pnorm(rnorm(100)), runif(100))    


## @knitr fig.keep="none"
qqplot(pt(rt(100,df=5),df=5),runif(100))    


## @knitr eval=FALSE
## shapiro.test(c(rnorm( 100), 5))
## shapiro.test(c(rnorm(1000), 5))
## shapiro.test(c(rnorm(4000), 5))


## @knitr 
f <- function(n,outlier=5) shapiro.test(c(rnorm(n),outlier))$p.value
M <- 500
res.100  = replicate(M, f(n=100))
res.1000 = replicate(M, f(n=1000))
res.4000 = replicate(M, f(n=4000))
sum(res.100  <= 0.05)/M
sum(res.1000 <= 0.05)/M
sum(res.4000 <= 0.05)/M


## @knitr 
res.100.nooutlier <- replicate(M, f(n=100, outlier=rnorm(1)))
sum(res.100.nooutlier <= 0.05)/500


