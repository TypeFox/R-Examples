
xshewhartrunsrules.ad <- function(mu1, mu0=0, c=1, type="12") {

# Shewhart chart
  if (type=="1") {
    p0 <- pnorm(  3*c, mean=mu1 ) - pnorm( -3*c, mean=mu1)
    ad <- 1/(1-p0)
  }

# ditto with runs rules
  if (type!="1") {
    Q1 <- xshewhartrunsrules.matrix(mu1, c=c, type=type)
    dimQ <- nrow(Q1)
    one <- rep(1, dimQ)
    I   <- diag(1, dimQ)
    arls <- solve(I-Q1, one)

    Q0 <- xshewhartrunsrules.matrix(mu0, c=c, type=type)
    psi <- Re(eigen(t(Q0))$vectors[,1])

    ad <- sum(psi * arls)/sum(psi)
  }

  ad
}