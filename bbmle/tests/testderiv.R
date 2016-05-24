library(bbmle)
old_opt <- options(digits=3)
## source("../R/dists.R")
## source("../R/mle.R")

## an attempt to sketch out by hand
##  how one would derive an analytic
##  gradient function for a formula-specified
##  likelihood and use it ...

## chain rule should be:

## deriv(probability distribution)/[prob params] *
##    deriv([prob params])/[model params] *
##   {OPTIONAL} deriv([model params])/[linear model params]

set.seed(1001)
x <- rbinom(50,size=10,prob=0.4)
suppressWarnings(mle2(x~dbinom(prob=p,size=10),start=list(p=0.3),data=data.frame(x)))

## step 1: construct gradient function for simplest example
f <- sbinom(prob=0.1,size=1)$formula

d1 <- deriv(parse(text=f),"prob",function.arg=TRUE)

## step 2: chain rule step #1
mle2(x~dbinom(prob=plogis(logitp),size=10),start=list(logitp=-1),
     data=data.frame(x))

f <- sbinom(prob=NA,size=NA)$formula

## note: plogis is not in derivatives table!!
##  will need to extend by text substitution ...
gsub("plogis(\\([^)]+\\))",
     "(1+exp(\\1))^(-1)",
     "plogis(logitprob)")

f2 <- gsub("plogis(\\([^)]+\\))",
     "(1+exp(\\1))^(-1)","plogis(logitp)")

## start with a single parameter (ignore 'size')
fun1 <- deriv(parse(text=f),c("prob"),function.arg=TRUE)
fun2 <- deriv(parse(text=f2),"logitp", function.arg=TRUE)

size <- 10
a1 <- attr(fun2(logitp=0),"gradient")
a2 <- attr(fun1(prob=plogis(0)),"gradient")

## compute gradient by variable and sum
colSums(apply(a1,2,"*",a2))
## rep(a1,length(x))*a2


## eventually we will want to do something tricky to
## 'memoise' results because optim() requires the
## objective function and gradient to be computed
## *separately*.  Not worth worrying about this in the
## first pass!
options(old_opt)
