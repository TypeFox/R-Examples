library(bbmle)

## Simulate data

set.seed(1)
x <- 1:5
y <- 2*x+1
noise <- rnorm(5, 0, 0.1)
mydata <- data.frame(x = x, y=y+noise)

## Model definition

model <- function(a, b) with(mydata, a*x+b)

## Negative log-likelihood

nll <- function(par) with(mydata, {
  a <- par[1]
  b <- par[2]
  sum(0.5*((y-model(a,b))/0.1)^2)
  
})

gr <- function(par) with(mydata, {
  a <- par[1]
  b <- par[2]
  dnllda <- -sum(((y-model(a,b))/0.1)*x/0.1)
  dnlldb <- -sum(((y-model(a,b))/0.1)*1/0.1)
  return(c(dnllda, dnlldb))
  
})

## optimization

parnames(nll) <- c("a", "b")
parnames(gr) <- c("a", "b")

fit <- mle2(nll, c(a = 1, b=2), gr=gr)

myprof <- profile(fit)

fit <- mle2(nll, c(a = 1, b=2), gr=gr, skip.hessian=TRUE)
myprof2 <- profile(fit,std.err=c(0.1,0.1))
