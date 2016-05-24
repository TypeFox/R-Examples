prop( ~ (diffmean <= -0.79 | diffmean >= 0.79), data = Randomization.Smiles)
2 * prop( ~ diffmean >= 0.79, data = Randomization.Smiles )
dotPlot( ~ diffmean, width = .03, cex = .5, groups = (diffmean >= 0.79), data = Randomization.Smiles)

