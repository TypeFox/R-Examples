Rand.toys <- do(5000) * rflip(16, .5)
head(Rand.toys, 3)
2 * prop( ~ (prop >= p.hat), data = Rand.toys)
histogram(~prop, data = Rand.toys, width=1/16)

