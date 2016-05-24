anova(Ants.Model)
MSE<- 138.7
mean(Ants ~ Filling, data = SandwichAnts)
mean <- 34.0
t.star <- qt(.975, df = 21); t.star
mean - t.star * (sqrt(MSE) / sqrt(8))
mean + t.star * (sqrt(MSE) / sqrt(8))

