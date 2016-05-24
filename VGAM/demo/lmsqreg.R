# Demo for lmsqreg
# At the moment this is copied from lms.bcn.Rd


data(bmi.nz, package = "VGAM")
fit <- vgam(BMI ~ s(age, df = c(4, 2)), lms.bcn(zero = 1),
            data = bmi.nz, trace = TRUE)
head(predict(fit), 3)
head(fitted(fit), 3)
head(bmi.nz, 3)
# Person 1 is near the lower quartile of BMI amongst people his age
head(cdf(fit), 3)

# Quantile plot
par(bty = "l", mar = c(5, 4, 4, 3) + 0.1, xpd = TRUE, mfrow = c(1, 2))
qtplot(fit, percentiles = c(5, 50, 90, 99), main = "Quantiles",
       xlim = c(15, 90), las = 1, ylab = "BMI", lwd = 2, lcol = 4)

# Density plot
ygrid <- seq(15, 43, len = 100)  # BMI ranges
par(lwd = 2)
aa <- deplot(fit, x0 = 20, y = ygrid,
           main = "Density functions at Age = 20, 42 and 55", xlab = "BMI")
aa
aa <- deplot(fit, x0 = 42, y = ygrid, add = TRUE, lty = 2, col = "orange")
aa <- deplot(fit, x0 = 55, y = ygrid, add = TRUE, lty = 4, col = 4, Attach = TRUE)
aa@post$deplot  # Contains density function values


