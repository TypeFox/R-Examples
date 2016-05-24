chisq.sample <- do(1000) *  chisq( ~ resample(toupper(letters[1:5]), 400))
histogram(~X.squared, data = chisq.sample)

