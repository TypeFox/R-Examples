# nbr2_6_2_3.r  Chanages to multivariable models
# Table 6.7  Hilbe,JM (2011), Negative Binomial Regression, 2 ed, Cambridge Univ Press
nobs <- 50000
x1m <- runif(nobs)
ym <-rpois(nobs, exp(1 + 0.5*x1 -.25*x2))
poim <- glm(ym ~ x1 + x2, family=poisson)
summary(poim)
mum <- predict(poim, type="response")
# change response ym
ym10 <- ym*10
poim10 <- glm(ym10 ~ x1 + x2, family=poisson)
summary(poim10)
mum10 <- predict(poim10, type="response")
# change predictor x1
xm10 <- x1*10
poix10 <- glm(ym ~ xm10 + x2, family=poisson)
summary(poix10)
mumx10 <- predict(poix10, type="response")



