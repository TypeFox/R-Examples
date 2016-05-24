Rand.Ants <- do(1000) * anova(lm(Ants ~ shuffle(Filling), data = SandwichAnts))
tally( ~ (F >= 5.63), data = Rand.Ants)
prop( ~ (F >= 5.63), data = Rand.Ants)
dotPlot( ~ F, width = 0.20, groups = (F <=5.63), data = Rand.Ants)

