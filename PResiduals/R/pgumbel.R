pGumbel <- function (q, loc = 0, scale = 1, lower.tail = TRUE) 
{
    q <- (q - loc)/scale
    p <- exp(-exp(q))
    if (lower.tail) 
        1 - p
    else p
}

pgumbel <- function(q, loc = 0, scale = 1, lower.tail = TRUE)
{
    q <- (q - loc)/scale
    p <- exp(-exp(-q))
    if (!lower.tail) 1 - p else p
}

dgumbel <- function (x, loc = 0, scale = 1, log = FALSE)
{
    x <- (x - loc)/scale
    d <- log(1/scale) - x - exp(-x)
    if (!log) ifelse(is.finite(x), exp(d), 0) else ifelse(is.finite(x), d, -Inf)
}

qgumbel <- function(p, loc = 0, scale = 1, lower.tail = TRUE)
{
  if(!lower.tail) p <- 1 - p 
  q <- -log(-log(p))
  q <- loc + scale * q
  return(q)
}

