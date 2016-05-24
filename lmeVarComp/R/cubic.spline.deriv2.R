cubic.spline.deriv2 <-
function(x, knots, interval)
{
  bspline.deriv(x, 4L, knots, interval, 2L)
}
