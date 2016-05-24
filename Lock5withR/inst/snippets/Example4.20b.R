# a 2-sided p-value is the sum of the values above
prop(~ (prop <= 0.2 | prop >= 0.8), data = RandomizationDist)
# We can also approximate the p-value by doubling one side
2 * prop(~ prop >= 0.80, data = RandomizationDist)

