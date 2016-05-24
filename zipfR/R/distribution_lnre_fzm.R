tdlnre.lnre.fzm <- function (model, x, ...)
{
  if (! inherits(model, "lnre.fzm")) stop("first argument must be object of class 'lnre.fzm'")
  
  alpha <- model$param$alpha
  A <- model$param$A
  B <- model$param$B
  C <- model$param2$C
  
  below <- x < A
  above <- x > B
  ok <- !(below | above)

  d <- numeric(length(x))
  d[below] <- 0
  d[above] <- 0
  d[ok] <- C * (x[ok] ^ (-alpha - 1))

  d
}

  
tplnre.lnre.fzm <- function (model, q, lower.tail=FALSE, ...)
{
  if (! inherits(model, "lnre.fzm")) stop("first argument must be object of class 'lnre.fzm'")
  
  alpha <- model$param$alpha
  A <- model$param$A
  B <- model$param$B
  C <- model$param2$C
  S <- model$S
  
  below <- q < A
  above <- q > B
  ok <- !(below | above)

  p <- numeric(length(q))
  p[below] <- S
  p[above] <- 0
  p[ok] <- (C / alpha) * (q[ok]^(-alpha) - B^(-alpha))

  if (lower.tail) p <- S - p
  p
}

  
tqlnre.lnre.fzm <- function (model, p, lower.tail=FALSE, ...)  
{
  if (! inherits(model, "lnre.fzm")) stop("first argument must be object of class 'lnre.fzm'")
  
  alpha <- model$param$alpha
  A <- model$param$A
  B <- model$param$B
  C <- model$param2$C
  S <- model$S

  if (lower.tail) p <- S - p
  
  below <- p < 0
  above <- p > S
  ok <- !(below)

  q <- numeric(length(p))
  q[below] <- B
  q[above] <- A
  q[ok] <- ( (alpha / C) * p[ok] + B^(-alpha) ) ^ (-1 / alpha)

  q
}


dlnre.lnre.fzm <- function (model, x, ...)
{
  if (! inherits(model, "lnre.fzm")) stop("first argument must be object of class 'lnre.fzm'")
  
  alpha <- model$param$alpha
  A <- model$param$A
  B <- model$param$B
  C <- model$param2$C
  
  below <- x < A
  above <- x > B
  ok <- !(below | above)

  d <- numeric(length(x))
  d[below] <- 0
  d[above] <- 0
  d[ok] <- C * x[ok]^(-alpha)

  d
}


plnre.lnre.fzm <- function (model, q, lower.tail=TRUE, ...)
{
  if (! inherits(model, "lnre.fzm")) stop("first argument must be object of class 'lnre.fzm'")
  
  alpha <- model$param$alpha
  A <- model$param$A
  B <- model$param$B
  C <- model$param2$C
  
  below <- q < A
  above <- q > B
  ok <- !(below | above)

  p <- numeric(length(q))
  p[below] <- 0
  p[above] <- 1
  p[ok] <- (q[ok]^(1-alpha) - A^(1-alpha)) / (B^(1-alpha) - A^(1-alpha))

  if (!lower.tail) p <- 1 - p
  p
}


qlnre.lnre.fzm <- function (model, p, lower.tail=TRUE, ...)
{
  if (! inherits(model, "lnre.fzm")) stop("first argument must be object of class 'lnre.fzm'")
  if (lower.tail) p <- 1 - p
  
  alpha <- model$param$alpha
  A <- model$param$A
  B <- model$param$B
  C <- model$param2$C

  below <- p < 0
  above <- p > 1
  ok <- !(below | above)

  q <- numeric(length(p))
  q[below] <- A
  q[above] <- B
  q[ok] <- ( ((1-alpha)/C) * p[ok] + A^(1-alpha) ) ^ (1/(1-alpha))

  q
}
