# Computes BIC for the regression of y on x, after doing subset selection 
# on the predictors
# require(BMA)

BICreg <- function(x, y)
{
  x <- as.matrix(x)
  y <- as.vector(y)
  n <- length(y)
  mod <- bicreg(y = y, x = x, nbest = 1)
  subset <- which(mod$which[1,])
  mod <- lm.fit(y = y, x = cbind(1,x[,subset]))
  # calculate the BIC for the regression
  sigma <- sqrt(mean(residuals(mod)^2))
  p <- n - df.residual(mod) + 1
  -n*log(2*pi) -2*n*log(sigma) -n -log(n)*p
}