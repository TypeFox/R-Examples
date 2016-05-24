# Demo for vgam


data(hunua, package = "VGAM")
fit.h <- vgam(agaaus ~ s(altitude), binomialff, data = hunua)
plot(fit.h, se = TRUE, lcol = "blue", scol = "orange", llwd = 2,
     slwd = 2, las = 1)


nn <- nrow(hunua)
ooo <- with(hunua, order(altitude))
with(hunua, plot(altitude[ooo], fitted(fit.h)[ooo], type = "l",
                 ylim = 0:1, lwd = 2, col = "blue", las = 1))
points(agaaus + (runif(nn)-0.5)/30 ~ altitude, hunua, col = "orange")


