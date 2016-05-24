tdlnre.lnre.zm <- function (model, x, ...)
{
  if (! inherits(model, "lnre.zm")) stop("first argument must be object of class 'lnre.zm'")
  
  alpha <- model$param$alpha
  B <- model$param$B
  C <- model$param2$C
  
  below <- x < 0
  above <- x > B
  ok <- !(below | above)

  d <- numeric(length(x))
  d[below] <- 0
  d[above] <- 0
  d[ok] <- C * (x[ok] ^ (-alpha - 1))

  d
}

  
tplnre.lnre.zm <- function (model, q, lower.tail=FALSE, ...)
{
  if (! inherits(model, "lnre.zm")) stop("first argument must be object of class 'lnre.zm'")
  if (lower.tail)
    stop("lower-tail type distribution not meaningful for LNRE model with S = Inf")
  
  alpha <- model$param$alpha
  B <- model$param$B
  C <- model$param2$C
  S <- Inf
  
  below <- q < 0
  above <- q > B
  ok <- !(below | above)

  p <- numeric(length(q))
  p[below] <- Inf
  p[above] <- 0
  p[ok] <- (C / alpha) * q[ok]^(-alpha) - (1 - alpha) / (B * alpha)

  p
}

  
tqlnre.lnre.zm <- function (model, p, lower.tail=FALSE, ...)  
{
  if (! inherits(model, "lnre.zm")) stop("first argument must be object of class 'lnre.zm'")
  if (lower.tail)
    stop("lower-tail type distribution not meaningful for LNRE model with S = Inf")
  
  alpha <- model$param$alpha
  B <- model$param$B
  C <- model$param2$C

  below <- p < 0
  ok <- !(below)

  q <- numeric(length(p))
  q[below] <- B
  q[ok] <- ( (alpha / C) * p[ok] + B^(-alpha) ) ^ (-1 / alpha)

  q
}


dlnre.lnre.zm <- function (model, x, ...)
{
  if (! inherits(model, "lnre.zm")) stop("first argument must be object of class 'lnre.zm'")
  
  alpha <- model$param$alpha
  B <- model$param$B
  C <- model$param2$C
  
  below <- x < 0
  above <- x > B
  ok <- !(below | above)

  d <- numeric(length(x))
  d[below] <- 0
  d[above] <- 0
  d[ok] <- C * x[ok]^(-alpha)

  d
}


plnre.lnre.zm <- function (model, q, lower.tail=TRUE, ...)
{
  if (! inherits(model, "lnre.zm")) stop("first argument must be object of class 'lnre.zm'")
  
  alpha <- model$param$alpha
  B <- model$param$B
  C <- model$param2$C
  
  below <- q < 0
  above <- q > B
  ok <- !(below | above)

  p <- numeric(length(q))
  p[below] <- 0
  p[above] <- 1
  p[ok] <- (q[ok] ^ (1-alpha)) / (B ^ (1-alpha))

  if (!lower.tail) p <- 1 - p
  p
}


qlnre.lnre.zm <- function (model, p, lower.tail=TRUE, ...)
{
  if (! inherits(model, "lnre.zm")) stop("first argument must be object of class 'lnre.zm'")
  if (lower.tail) p <- 1 - p
  
  alpha <- model$param$alpha
  B <- model$param$B
  C <- model$param2$C

  below <- p < 0
  above <- p > 1
  ok <- !(below | above)

  q <- numeric(length(p))
  q[below] <- 0
  q[above] <- B
  q[ok] <- B * (p[ok] ^ (1 / (1-alpha)))

  q
}
