`quadsol` <-
function(a, b, c)
{
  if(a > 0) {
    # Return the smallest nonnegative integer n such that an^2+bn+c>=0
    delta <- b * b - 4 * a * c
    if(delta < 0) {
      0
    }
    else {
      root1 <- (-b-sqrt(delta))/2/a
      root2 <- (-b+sqrt(delta))/2/a
      ifelse(root1 >= 0, 0, max(0, ceiling(root2)))
    }
  }
}

