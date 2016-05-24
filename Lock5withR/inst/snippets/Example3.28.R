median(~ Price, data = MustangPrice )
Boot.Mustang <- do(5000) * median( ~Price, data = resample(MustangPrice) )
head(Boot.Mustang, 3)
histogram( ~ median, n = 50, data = Boot.Mustang)

