RandomizationDist <- do(1000) * rflip(10, .5)  # 10 because n = 10
head(RandomizationDist)
histogram( ~ prop, label = TRUE, width = 1/10, data = RandomizationDist)

