partialCorrelation <- function(x) {
  conc <- solve(var(x))
  resid.sd <- 1/sqrt(diag(conc))
  - sweep(sweep(conc, 1, resid.sd, "*"), 2, resid.sd, "*")
}

## cor(x,y)-cor(x,z)*cor(y,z))/(sqrt(1-cor(x,z)^2)*sqrt(1-cor(y,z)^2))

