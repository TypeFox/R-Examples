

# ========================================================================
# rvhyper  -  hypergeometric rvs
# ========================================================================

rvhyper <- function(nn=1, m, n, k) {
  warning("NOT YET READY")
  attr(nn, "n.name") <- "nn"
  rvvapply(stats:::rhyper, n.=nn, m=m, n=n, k=k)
}

