BootP <- do(1000) * rflip(100, .52)
head(BootP, 3)
dotPlot(~ prop, width = .01, data = BootP)

