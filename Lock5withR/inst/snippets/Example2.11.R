ICU20 <- subset(ICUAdmissions, Age == "20")
mean( ~ HeartRate, data = ICU20)
median( ~ HeartRate, data = ICU20)
ICU55 = subset(ICUAdmissions, Age == "55")
mean( ~ HeartRate, data = ICU55)
median( ~ HeartRate, data = ICU55)

