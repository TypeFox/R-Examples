## ---- echo=FALSE,message=FALSE, warning=FALSE----------------------------
#devtools::load_all()

## ----fit-----------------------------------------------------------------
# middle
library(cdfquantreg)
data(cdfqrExampleData)
ipcc_mid <- subset(IPCC, mid == 1 & high == 0)

# We use T2-T2 distribution
fd <- "t2"
sd <- "t2"

# Fit the null model
fit_null <- cdfquantreg(probm ~ 1 | 1, fd, sd, data = ipcc_mid)

# Fit the target model
fit <- cdfquantreg(probm ~ valence | valence, fd, sd, data = ipcc_mid)

# Obtain the statistics for the null model
summary(fit)

## ----plotfit-------------------------------------------------------------
# Compare the empirical distribution and the fitted values distribution
plot(fit)

# Plot the fitted values
plot(fitted(fit, "full"))

# Check Residuals
plot(residuals(fit, "raw"))

