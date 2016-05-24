Boot.Temp <- do (5000) * mean(~ BodyTemp, data = resample(BodyTemp50))
head(Boot.Temp,3)
mean(~ mean, data = Boot.Temp)
cdata( ~ mean, 0.95, data = Boot.Temp)
histogram( ~ mean, width = .01, v = c(98.26, 98.6), 
          groups = (98.05 <= mean & mean <= 98.46), data = Boot.Temp)

