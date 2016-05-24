# Now we'll do it 1000 times
Bootstrap <- do(1000) * mean( ~Time, data = resample(CommuteAtlanta))
head(Bootstrap, 3)
# We should check that that our bootstrap distribution has an appropriate shape:
dotPlot( ~ mean, width = 0.005, data = Bootstrap)

