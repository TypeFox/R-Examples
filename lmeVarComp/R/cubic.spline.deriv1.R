cubic.spline.deriv1 <-
function(x, knots, interval)
{
  bspline.deriv(x, 4L, knots, interval, 1L)
}
