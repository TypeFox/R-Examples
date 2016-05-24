BootM <- do(5000) * cor(Price~Miles, data = resample((MustangPrice)))
head(BootM, 3)

