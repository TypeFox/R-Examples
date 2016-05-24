Randomization.Temp2 <- 
  do(5000) * ( mean( ~ BodyTemp, data = resample(BodyTemp50) ) - 98.26 ) 
head(Randomization.Temp2, 3)
mean( ~ result, data = Randomization.Temp2 )

