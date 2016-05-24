# nbr2_6_2_2.r  Chanages to univariable models
# Table 6.6  Hilbe,JM (2011), Negative Binomial Regression, 2 ed, Cambridge Univ Press
nobs <- 50000
x1 <- runif(nobs)
y <-rpois(nobs, exp(1 + 0.5*x1))
poi <- glm(y ~ x1, family=poisson)
summary(poi)
mu <- predict(poi, type="response")
# change response y
y10 <- y*10
poi10 <- glm(y10 ~ x1, family=poisson)
summary(poi10)
mu10 <- predict(poi10, type="response")
# change predictor x1
x10 <- x1*10
poix <- glm(y ~ x10, family=poisson)
summary(poix)
mux10 <- predict(poix, type="response")
