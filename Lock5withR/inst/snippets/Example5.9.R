Randomization.Temp <- do(10000) * (mean( ~ BodyTemp, data = resample(BodyTemp50)) + 0.34)
histogram( ~ mean, width = .025, fit = "normal", data = Randomization.Temp)

