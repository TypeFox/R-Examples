"ellipse.profile.nls" <-
  function (x, which = c(1, 2), level = 0.95, t = sqrt(2 * qf(level, 2, attr(x, 
                               "summary")$df[2])), npoints = 100, ...) 
{
  ellipse.profile(x, which = which, level = level, t = t, npoints = npoints)
}
