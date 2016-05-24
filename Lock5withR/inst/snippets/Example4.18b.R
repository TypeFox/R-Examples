Randomization.Caff <- do (1000) * ediff( mean( Taps ~ shuffle(Group), data = CaffeineTaps ) )
head(Randomization.Caff,3)
dotPlot( ~ No.Caffeine, width = 0.2, data = Randomization.Caff)

