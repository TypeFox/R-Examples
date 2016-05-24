descStat <- function (x, na.rm = TRUE)
{
  if(!is.numeric(x))
    stop("'x' must be numeric")
  m <- mean(x, na.rm = na.rm)
  s <- sd(x, na.rm = na.rm)
  c(n = length(x),
    obs = sum(!is.na(x)),
    mean = m,
    median = median(x, na.rm = na.rm),
    sd = s,
    cv = s/m,
    min = min(x, na.rm = na.rm),
    max = max(x, na.rm = na.rm))
}

