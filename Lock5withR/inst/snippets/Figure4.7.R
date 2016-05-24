prop( ~ (prop >= 0.64), data = Randomization.Match ) # 16/25
histogram( ~ prop, width = 0.04, groups = (prop >= 0.64), data = Randomization.Match)

