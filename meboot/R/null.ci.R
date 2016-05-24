
null.ci <- function(x, level = 0.95, null.value = 0, type = 8, ...) 
{
  q <- mean(x <= null.value) * (1 - level)
  quantile(x, c(q, q + level), type = type, ...)
}
