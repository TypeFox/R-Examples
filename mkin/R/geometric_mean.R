geometric_mean <- function(x, na.rm = FALSE)
{
  exp(mean(log(x), na.rm = na.rm))
}
