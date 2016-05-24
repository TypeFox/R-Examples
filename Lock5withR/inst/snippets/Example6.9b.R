Randomization <- do(1000) * rflip(9, 0.5)
head(Randomization, 3)
prop( ~ (prop >= p.hat), data = Randomization)

