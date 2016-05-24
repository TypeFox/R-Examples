data(sp.raw)
data(spto87)

tsp.raw <- attributes(sp.raw)$times
tspto87 <- attributes(spto87)$times
plot(tsp.raw, sp.raw, type = "l")
plot(tspto87, -spto87, type = "l")

out1 <- gev(-spto87, "year")
out2 <- gev(-spto87, "semester")
plot(tspto87, -spto87, type = "l", ylim = c(-6, 15))
rlevel.gev(out1, 10, add = TRUE)
rlevel.gev(out2, 20, add = TRUE)

plot(tspto87, -spto87, type = "l", ylim = c(-6, 25))
rlevel.gev(out1, 50, add = TRUE)
n <- length(sp.raw)
sp <- log(sp.raw[-1]/sp.raw[-n])
plot(tsp.raw[-1], -sp, type = "l")
