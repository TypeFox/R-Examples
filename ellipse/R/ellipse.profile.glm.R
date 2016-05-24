"ellipse.profile.glm" <-
  function (x, which = c(1, 2), level = 0.95, t, npoints = 100, dispersion, ...) 
{
  if (missing(dispersion)) dispersion <- ifelse(attr(x,"original.fit")$family$family %in% c('poisson','binomial'),
			                                      1, NA)
  if (missing(t)) t <- ifelse(is.na(dispersion),sqrt(2 * qf(level, 2, attr(x,"summary")$df[2])),
			                                    sqrt(qchisq(level, 2)*dispersion/attr(x,"summary")$dispersion))
  ellipse.profile(x, which = which, level = level, t = t, npoints = npoints)
}
