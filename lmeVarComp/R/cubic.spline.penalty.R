cubic.spline.penalty <-
function(knots, interval, derivative = 2L)
{
  D <- if (derivative == 1L) {
    cubic.spline.penalty1(knots, interval)
  } else {
    cubic.spline.penalty2(knots, interval)
  }
  
  D
}
