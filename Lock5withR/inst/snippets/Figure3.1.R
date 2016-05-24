# Now we'll do it 1000 times
Sampledist <- do(1000) * mean( ~ FTGradEnrollment, data = sample(StatisticsPhD, 10))
head(Sampledist, 3)
dotPlot( ~ mean, width = .005, data = Sampledist)

