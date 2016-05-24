fit <- lm(metabolic.rate ~ body.weight, data=rmr)
summary(fit)
811.2267 + 7.0595 * 70 # , or:
predict(fit, newdata=data.frame(body.weight=70))
qt(.975,42)
7.0595 + c(-1,1) * 2.018 * 0.9776 # , or:
confint(fit)
summary(lm(log(ab)~age, data=malaria))
plot(log(ab)~age, data=malaria)
rho <- .90 ; n <- 100
x <- rnorm(n)
y <- rnorm(n, rho * x, sqrt(1 - rho^2))
plot(x, y)
cor.test(x, y)
cor.test(x, y, method="spearman")
cor.test(x, y, method="kendall")
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
