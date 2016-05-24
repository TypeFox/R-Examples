treb.model <- lm(distance~projectileWt,data=trebuchet2)
ef <- treb.model$effects; n <- length(ef)
ef
# total length
sum(ef^2)
sum(trebuchet2$distance^2)
# projection of residuals into n-2 orthogonal components
sum(ef[3:n]^2)
sum(resid(treb.model)^2)
# projection in direction of u[0] is mean * sqrt(n)
# Note: R's u[0] points in the opposite direction.
mean(trebuchet2$distance) * sqrt(n)
ef[1]
# beta1.hat obtained from projection in direction of u[1]
v1 <- trebuchet2$projectileWt - mean(trebuchet2$projectileWt)
ef[2] / sqrt(sum(v1^2))
coef(treb.model)
