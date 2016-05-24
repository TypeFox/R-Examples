# Demo for binom2.or



data(hunua, package = "VGAM")
Hunua <- hunua
Hunua <- transform(Hunua, y00 = (1-agaaus) * (1-kniexc),
                          y01 = (1-agaaus) *    kniexc,
                          y10 =    agaaus  * (1-kniexc),
                          y11 =    agaaus  *    kniexc)



fit <- vgam(cbind(y00, y01, y10, y11) ~ s(altitude, df = c(4, 4, 2.5)),
            binom2.or(zero = NULL), data = Hunua)
par(mfrow = c(2, 3))
plot(fit, se = TRUE, scol = "darkgreen", lcol = "blue")
summary(fit)


# Plot the marginal functions together
mycols <- c("blue", "orange")
plot(fit, which.cf = 1:2, lcol = mycols, scol = mycols,
     overlay = TRUE, se = TRUE, llwd = 2, slwd = 2)
legend(x = 100, y = -4, leg = c("Agathis australis", "Knightia excelsa"),
       col = mycols, lty = 1)


# Plot the odds ratio
ooo <- order(fit@x[, 2])
plot(fit@x[ooo, 2], exp(predict(fit)[ooo, "log(oratio)"]), 
     log = "y", xlab = "Altitude (m)", ylab = "Odds ratio (log scale)",
     col = "blue", type = "b", las = 1)
abline(h = 1, lty = 2)  # Denotes independence between species



