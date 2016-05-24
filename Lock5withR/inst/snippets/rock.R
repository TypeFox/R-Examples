p.hat <- 66/119; p.hat
p.hat * 119               # check >= 10
(1 - p.hat) * 119         # check >= 10
SE <- sqrt(1/3 * 2/3 / 119); SE
z <- (p.hat - 1/3) / SE; z
pnorm(z)                  # large side (rounded)
1 - pnorm(z)              # small side (less rounding)
2 * (1 - pnorm(z))        # p-value = 2 * small side

