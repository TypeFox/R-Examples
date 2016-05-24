Boot.Rent <- do(1000) * mean( ~ Rent, data = resample(ManhattanApartments)) 
head(Boot.Rent, 3)
favstats( ~ mean, data = Boot.Rent)
cdata( ~ mean, 0.95, data = Boot.Rent)

