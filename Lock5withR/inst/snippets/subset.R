head(ICUAdmissions, 2)
tally( ~ sex, data = ICUAdmissions) 
ICUMales <- subset(ICUAdmissions, sex == "Male") # notice the double =
tally( ~ sex, data = ICUMales)

