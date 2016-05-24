prop( ~ abs(result) > 0.34, data = Randomization.Temp2 )
histogram( ~ result, width = .01, v = c(0.34, -0.34), data = Randomization.Temp2)

