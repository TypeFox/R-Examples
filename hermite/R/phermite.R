 phermite <- function(q, a, b, m=2, lower.tail=TRUE)
{
  if (a > 20 | b > 20)     res <- ifelse(q>=0, edg(floor(q), a, b, m), 0)
  if (a <= 20 & b <= 20)   res <- ifelse(q>=0, vapply(q, function(x) sum(dhermite(0:x, a, b, m)), 1), 0)
  if (lower.tail == FALSE) res <- 1 - res
  return(res)
}
