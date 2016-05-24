prop( ~ (mean <= 30), data = Bootstrap)                # proportion less than 30 min
prop( ~ (mean >= 31), data = Bootstrap)                # proportion greater than 31 min
prop( ~ (mean >= 30 & mean <= 31), data = Bootstrap)   # proportion between 30 and 31 min

