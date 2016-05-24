z30 <- (30 - 29.11) / 0.93; z30  # z-score for 30 min
z31 <- (31 - 29.11) / 0.93; z31  # z-score for 31 min
xpnorm(c(30, 31), 29.11, 0.93)   # original normal distribution proportion between 30 and 31 min
xpnorm(c(z30, z31))              # standardized distribution proportion between 30 and 31 min
pnorm(z31) - pnorm(z30)

