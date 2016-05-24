
xshewhartrunsrules.arl <- function(mu, c=1, type="12") {

# Shewhart chart
  if (type=="1") {
    p0  <- pnorm(  3*c, mean=mu ) - pnorm( -3*c, mean=mu)
    arls <- 1/(1-p0)
  }

# ditto with runs rules
  if (type!="1") {
    Q <- xshewhartrunsrules.matrix(mu, c=c, type=type)
    dimQ <- nrow(Q)
    one <- rep(1, dimQ)
    I   <- diag(1, dimQ)
    arls <- solve(I-Q, one)
  }

  arl <- arls[1]
  arl
}