ltdlnre.lnre <- function (model, x, base=10, log.x=FALSE, ...)
{
  if (! inherits(model, "lnre")) stop("first argument must be object of class 'lnre'")

  if (log.x) x <- base ^ x
  log(base) * x * tdlnre(model, x)
}


ldlnre.lnre <- function (model, x, base=10, log.x=FALSE, ...)
{
  if (! inherits(model, "lnre")) stop("first argument must be object of class 'lnre'")

  if (log.x) x <- base ^ x
  log(base) * x * dlnre(model, x)
}


rlnre.lnre <- function (model, n, ...)
{
  if (! inherits(model, "lnre")) stop("first argument must be object of class 'lnre'")
  if (length(n) > 1) n <- length(n)

  if (length(n) == 0 || n == 0) {
    factor()                            # return empty factor vector
  }
  else {
    x <- runif(n)
    y <- 1 + floor( tplnre(model, qlnre(model, x)) ) # 1 + floor( G( F^-1(x) ) )
    factor(y)
  }
}
