Randomization.Coc <- do (5000) * diff( prop( response ~ shuffle(treatment), data = Cocaine ) )
head(Randomization.Coc)

