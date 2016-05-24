"ellipse.lm" <-
  function (x, which = c(1, 2), level = 0.95, t = sqrt(2 * qf(level, 
                                                2, x$df.residual)), ...) 
{
  s <- summary(x)
  ellipse.default(s$sigma^2 * s$cov.unscaled[which, which], 
                  centre = x$coefficients[which], t = t, ...)
}
