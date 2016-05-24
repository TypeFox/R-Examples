##
## demo1sine.r - symbolic Regression on sine / cosine function
##

require(rgp)

data1 <- {
  x1 <- seq(0, 4*pi, length.out=201)
  x2 <- seq(0, 4*pi, length.out=201)
  y <- sin(x1) + cos(2*x2)
  data.frame(y=y, x1=x1, x2=x2)
}

mdl <- symbolicRegression(y ~ x1 + x2, data=data1,
                          functionSet=arithmeticFunctionSet,
                          stopCondition=makeTimeStopCondition(120))

## Plot prediction
plot(data1$y, col=1, type="l"); points(predict(mdl, newdata = data1), col=2, type="l")

## Extract best and worst individual:
bf <- mdl$population[[which.min(sapply(mdl$population, mdl$fitnessFunction))]]
wf <- mdl$population[[which.max(sapply(mdl$population, mdl$fitnessFunction))]]

