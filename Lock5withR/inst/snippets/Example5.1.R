Bootstrap <- do(1000) * mean( ~ Time, data = resample(CommuteAtlanta))
head(Bootstrap, 3)
histogram( ~ mean, density = TRUE, data = Bootstrap)
densityplot( ~ mean, data = Bootstrap)

