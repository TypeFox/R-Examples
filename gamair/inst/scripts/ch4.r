# R code for chapter 4 of Wood (2006) "GAMs: An Introduction with R"

## 4.11.1 Backfitting GAMs
## setup up some data...
n<-400;sig<-2
x0 <- runif(n, 0, 1);x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1);x3 <- runif(n, 0, 1)
f0 <- function(x) 2 * sin(pi * x)
f1 <- function(x) exp(2 * x)
f2 <- function(x) 0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10
f <- f0(x0) + f1(x1) + f2(x2)
e <- rnorm(n, 0, sig)
y <- f + e
x <- cbind(x0,x1,x2,x3)
edf <- c(3,2,7,1.01);m <- 4

## backfit...

f<-x*0;alpha<-mean(y);ok <- TRUE;rss0 <- 0
while (ok) { # backfitting loop
  for (i in 1:m) { # loop through the smooth terms
    ep <- y - rowSums(f[,-i]) - alpha
    b <- smooth.spline(x[,i],ep,df=edf[i])
    f[,i] <- predict(b,x[,i])$y
  }
  rss <- sum((y-rowSums(f))^2)
  if (abs(rss-rss0)<1e-6*rss) ok <- FALSE
  rss0 <- rss
}
