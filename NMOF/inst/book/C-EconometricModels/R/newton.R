# newton.R -- version 2010-12-29
cf <- c(5, 5, 5, 5, 5, 105) # cashflows
tm <- 1:6              # maturities
ytm <- 0.046           # the "true" yield
b0 <- sum(cf/((1 + ytm)^tm))
cf <- c(-b0, cf); tm <- c(0, tm)

r <- 0.1   # initial value for r
h <- 1e-8  # fin-diff step
dr <- 1    # change in r
while (abs(dr) > 1e-5) {
    g <- sum(cf/((1 + r)^tm))
    dg <- (sum(cf/((1 + r + h)^tm)) - g)/h
    dr <- g/dg
    print(r <- r - dr)
}