BootP400 <- do(1000) * rflip(400, .52)
head(BootP400, 3)
cdata( ~ prop, 0.95, data = BootP400)
dotPlot(~ prop, width = 0.005, groups = (0.472 <= prop & prop<= 0.568), data = BootP400)

