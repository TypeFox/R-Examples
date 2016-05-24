# Demo for Zero-Inflated Poisson

set.seed(111)
zdata <- data.frame(x2 = runif(nn <- 1000))
zdata <- transform(zdata, pstr01  = logit(-0.5 + 1*x2, inverse = TRUE),
                          pstr02  = logit( 0.5 - 1*x2, inverse = TRUE),
                          Ps01    = logit(-0.5       , inverse = TRUE),
                          Ps02    = logit( 0.5       , inverse = TRUE),
                          lambda1 =  loge(-0.5 + 2*x2, inverse = TRUE),
                          lambda2 =  loge( 0.5 + 2*x2, inverse = TRUE))
zdata <- transform(zdata, y1 = rzipois(nn, lambda = lambda1, pstr0 = Ps01),
                          y2 = rzipois(nn, lambda = lambda2, pstr0 = Ps02))

with(zdata, table(y1))  # Eyeball the data
with(zdata, table(y2))
with(zdata, stem(y2))
fit1 <- vglm(y1 ~ x2, zipoisson(zero = 1), data = zdata, crit = "coef")
fit2 <- vglm(y2 ~ x2, zipoisson(zero = 1), data = zdata, crit = "coef")
coef(fit1, matrix = TRUE)  # These should agree with the above values
coef(fit2, matrix = TRUE)  # These should agree with the above values


head(fit1@misc$pobs0)  # The estimate of P(Y=0)

coef(fit1)
coef(fit1, matrix=TRUE)
Coef(fit1)


