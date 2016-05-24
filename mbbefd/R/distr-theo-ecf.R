#new exposure curve function for some distributions
#limited expected value is implemented in actuar

ecunif <- function(x, min = 0, max =1)
{
  levunif(x, min = min, max = max) / munif(1, min = min, max = max)
}

ecbeta <- function(x, shape1, shape2)
{
  levbeta(x, shape1 = shape1, shape2 = shape2) / mbeta(1, shape1 = shape1, shape2 = shape2)
}
