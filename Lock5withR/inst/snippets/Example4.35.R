Randomization.Temp3 <- 
  do(5000) * ( mean( ~ BodyTemp, data = resample(BodyTemp50) ) + 0.14 ) 
head(Randomization.Temp3, 3)
mean( ~ result, data = Randomization.Temp3 )
cdata(~ result, 0.95, data = Randomization.Temp3)
histogram( ~ result, width = .01, v = c(98.26, 98.4), 
          groups = (98.19 <= result & result <= 98.62), 
          xlim = c(97.8, 99.0), data = Randomization.Temp3) # randomization
histogram( ~ mean, width = .01, v = c(98.26, 98.4), 
          groups = (98.05 <= mean & mean <= 98.46), 
          xlim = c(97.8, 99.0), data = Boot.Temp) # bootstrap

