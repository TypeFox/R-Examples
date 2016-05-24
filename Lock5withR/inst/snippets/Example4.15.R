prop( ~ (prop >= 0.60), data = Randomization.Match ) # 15/25
prop( ~ (prop >= 0.76), data = Randomization.Match ) # 19/25
histogram( ~ prop, width = 0.04, groups = (prop >= 0.60), data = Randomization.Match )
histogram( ~ prop, width = 0.04, groups = (prop >= 0.76), data = Randomization.Match )

