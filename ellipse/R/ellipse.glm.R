"ellipse.glm" <-
  function (x, which = c(1, 2), level = 0.95, t, npoints = 100, dispersion, ...) 
{
  s <- summary(x)
  est.disp <- missing(dispersion) & !(x$family$family %in% c('poisson','binomial'))
  if (missing(dispersion)) dispersion <- s$dispersion
  if (missing(t)) t <- ifelse(est.disp,sqrt(2 * qf(level, 2, s$df[2])),
			                           sqrt(qchisq(level, 2)))
  ellipse.default(dispersion * s$cov.unscaled[which, which], 
                  centre = x$coefficients[which], t = t, npoints = npoints, ...)
}
