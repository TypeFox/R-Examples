cubic.spline <-
function(x, knots, interval)
{
  bspline(x, 4L, knots, interval)
}
