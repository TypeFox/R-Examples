
xDshewhartrunsrulesFixedm.arl <- function(delta, c=1, m=100, type="12") {
  mus <- (1:m)*delta

# Shewhart chart
  if (type=="1") {
    p0  <- pnorm(  3*c, mean=mus ) - pnorm( -3*c, mean=mus)
    arls <- 1/(1-p0[m])
    for ( i in (m-1):1 ) arls <- 1 + p0[i]*arls
  }

# ditto with runs rules
  if (type!="1") {
    Q <- xshewhartrunsrules.matrix(mus[m], c=c, type=type)
    dimQ <- nrow(Q)
    one <- rep(1, dimQ)
    I   <- diag(1, dimQ)
    arls <- solve(I-Q, one)

    for ( i in (m-1):1 ) {
      Q <- xshewhartrunsrules.matrix(mus[i], c=c, type=type)
      arls <- 1 + (Q %*% arls)[,1]
    }
  }

  arl <- arls[1]
  arl
}