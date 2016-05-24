# Example 10.7

histochart(dbinom(0:41, 41, 0.5) ~ factor(0:41))

# Normal approximation
2 * pnorm(31, mean = 20.5, sd = 3.2, lower.tail = FALSE)

# Exact binomial test
binom.test(31, 41, p = 0.5)
