

rvcauchy <- function (n=1, location=0, scale=1) {
  rvvapply(stats:::rcauchy, n.=n, location=location, scale=scale)
}


