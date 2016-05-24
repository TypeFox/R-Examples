radialregression <- function(TE, seTE, k){
  
  x <-  1/seTE
  y <- TE/seTE

  reg <- lm(y ~ x)

  intercept <- coefficients(reg)[1]
  slope     <- coefficients(reg)[2]

  se.slope     <- 1/sd(sqrt(1/seTE^2))/sqrt(k-1)
  se.intercept <- se.slope*mean(1/seTE^2)

  res <- list(intercept=intercept,
              se.intercept=se.intercept,
              slope=slope,
              se.slope=se.slope,
              k=k,
              sigma=summary(reg)$sigma,
              r.squared=summary(reg)$r.squared)
  res
}
