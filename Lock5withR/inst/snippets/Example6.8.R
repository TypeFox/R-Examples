p.hat <- 66/119; p.hat
p <- 1/3; p
p * 119                   # check >= 10
(1 - p) * 119             # check >= 10
SE <- sqrt(p * (1 - p) / 119); SE
z <- (p.hat - p) / SE; z
pnorm(z)                  # large side (rounded)
1 - pnorm(z)              # small side (less rounding)
2 * (1 - pnorm(z))        # p-value = 2 * small side
prop.test(66, 119, p = 1/3)

