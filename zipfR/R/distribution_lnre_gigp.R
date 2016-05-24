tdlnre.lnre.gigp <- function (model, x, ...)
{
  if (! inherits(model, "lnre.gigp")) stop("first argument must be object of class 'lnre.gigp'")
  
  gamma <- model$param$gamma
  b <- model$param$B                    # original notation from Baayen (2001)
  c <- model$param$C

  C <- (2 / (b*c))^(gamma+1) / besselK(b, gamma+1)
  d <- C * x^(gamma-1) * exp(- x/c - (b*b*c)/(4*x))
  
  d
}


dlnre.lnre.gigp <- function (model, x, ...)
{
  if (! inherits(model, "lnre.gigp")) stop("first argument must be object of class 'lnre.gigp'")
  
  gamma <- model$param$gamma
  b <- model$param$B                    # original notation from Baayen (2001)
  c <- model$param$C

  C <- (2 / (b*c))^(gamma+1) / besselK(b, gamma+1)
  d <- C * x^gamma * exp(- x/c - (b*b*c)/(4*x))
  
  d
}
