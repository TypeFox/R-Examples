# Some testing for a random generator for log-concave density estimates
require("REBayes")
require("logcondens")
# Test problem for comparison of log concave estimators
n <- 200
x <- rgamma(n,5)
tf <- system.time(f <- medde(x, lambda = -0.5))
tg <- system.time(g <- logConDens(x,smoothed = FALSE))
plot(g)
lines(f)
y <- rmedde(200,f)
points(y,(y/y)*.05,col = 2, cex = 0.5)
