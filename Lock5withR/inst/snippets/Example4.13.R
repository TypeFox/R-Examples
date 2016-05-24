Randomization.Match <- 
  do(10000) * rflip(25, .5)  # 25 because n = 25
head(Randomization.Match)
histogram( ~ prop, width = 0.04, data = Randomization.Match)

