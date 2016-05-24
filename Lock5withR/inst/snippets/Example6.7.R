p.hat <- 0.19; p.hat
p <- 0.20; p
p * 1013                   # check >= 10
(1 - p) * 1013             # check >= 10
SE <- sqrt(p * (1 - p) / 1013); SE
z <- (p.hat - p) / SE; z
pnorm(z)

