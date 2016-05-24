## This illustration uses the ggplot2 package to construct plots that
## convey more information than the default plot function in base R

rm(list=ls())
library(crs)
set.seed(42)

## Simulate data then estimate a crs model having one continuous and
## one categorical predictor

n <- 1000
x <- runif(n)
neval <- 100
c <- 3
z <- as.integer(cut(runif(n),breaks=qunif(seq(0,1,length=c+1))))-1
dgp <- cos(4*pi*x)+z
dgp <- dgp/sd(dgp)
y <- dgp + rnorm(n,sd=.25)
z <- factor(z)

model <- crs(y ~ x + z)

## Create evaluation data 

data.eval <- expand.grid(x = seq(min(x), max(x), length = neval),
                    z = levels(z))

library(ggplot2)

## Create a plot using qplot() that will use colors for observations
## classified by levels of the factor and the same colors for the
## predicted model values

data.eval$y <- predict(model, newdata=data.eval)
qplot(x, y, colour=z) + geom_line(data=data.eval) 

## Create a plot using qplot() that will use colors for observations
## classified by levels of the factor and the same colors for the
## predicted model values but add

data.eval$ucl <- attr(data.eval$y,"upr")
data.eval$lcl <- attr(data.eval$y,"lwr")

qplot(x, y, colour=z) + geom_smooth(aes(ymin = lcl, ymax = ucl), data=data.eval, stat="identity") 
