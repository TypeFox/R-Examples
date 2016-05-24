data(euroex)
length(euroex)
plot(euroex, type = "l", xlab = "Day Number", ylab = "Exchange Rate")
euroex.ldr <- log(euroex[-1]/euroex[-975])
plot(euroex.ldr, type = "l", xlab = "Day Number", ylab = "Log-Daily Return")
euroex.ldr <- 100*euroex.ldr
mrl.plot(euroex.ldr)
sum(euroex.ldr > 0.9)
gpd.fitrange(euroex.ldr, -1, 1.4, nint = 100)
euroex.gpd <- gpd.fit(euroex.ldr, 0.9, npy = 250)
gpd.diag(euroex.gpd)
gpd.profxi(euroex.gpd, -0.5, 0.3)
gpd.prof(euroex.gpd, m = 10, npy = 250, 1.75, 3.5)
trend <- as.matrix((-500:473)/250)
euroex.gpd2 <- gpd.fit(euroex.ldr, 0.9, npy = 250, ydat = trend,
    sigl = 1, siglink = exp)
gpd.diag(euroex.gpd2)


