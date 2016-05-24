# Example 10.4

# Pr[Males under 148.6 cm or over 193 cm]
pnorm(148.6, mean = 175.6, sd = 7.1) + 
  pnorm(193, mean = 175.6, sd = 7.1, lower.tail = FALSE)

# Pr[Females under 148.6 cm or over 193 cm]
pnorm(148.6, mean = 162.6, sd = 6.4) + 
  pnorm(193, mean = 162.6, sd = 6.4, lower.tail = FALSE)
