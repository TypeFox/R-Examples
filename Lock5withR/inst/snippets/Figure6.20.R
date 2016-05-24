Boot.Love <- do(5000)*diff(prop(Response~Gender, data = resample(OneTrueLove)))
head(Boot.Love, 3)
histogram(~Agree.Male, fit = "normal", data = Boot.Love)

