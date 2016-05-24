prop( ~ abs(result - 98.6) > 0.34, data = Randomization.Temp)
histogram( ~ result, width = .01, v = c(98.40, 98.6, 98.81), data = Randomization.Temp)

