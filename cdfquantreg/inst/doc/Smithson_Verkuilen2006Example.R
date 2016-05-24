## ---- echo=FALSE, warning=FALSE, message=FALSE---------------------------
# library(devtools)
# load_all()

## ------------------------------------------------------------------------
library(cdfquantreg)
data(cdfqrExampleData)
#Quick overview of the variables
rbind(head(JurorData,4),tail(JurorData,4))

## ----fit-----------------------------------------------------------------
# We use T2-T2 distribution
fd <- "logit" # The parent distribution
sd <- "logistic" # The child distribution

# Fit the null model
fit_null <- cdfquantreg(crc99 ~ 1 | 1, fd, sd, data = JurorData)

# Fit a main effect model
fit1 <- cdfquantreg(crc99 ~ vert + confl | 1, fd, sd, data = JurorData)

# Fit the full model
fit2 <- cdfquantreg(crc99 ~ vert*confl | 1, fd, sd, data = JurorData)

anova(fit1,fit2)

# Obtain the statistics for the null model
summary(fit2)

## ----fit_dismod----------------------------------------------------------
# Fit a main effect model
fit3 <- cdfquantreg(crc99 ~ vert*confl |vert + confl, fd, sd, data = JurorData)

# Fit the full model
fit4 <- cdfquantreg(crc99 ~ vert*confl |vert*confl, fd, sd, data = JurorData)

anova(fit2, fit3, fit4)

# Obtain the statistics for the null model
summary(fit4)

## ----fig.height=4, fig.width=4-------------------------------------------
# Compare the empirical distribution and the fitted values distribution
breaks <- seq(0,1,length.out =11)

plot(fit4,xlim = c(0.1,1),ylim = c(0,3), breaks = breaks)

## ----plotfit,fig.height=9, fig.width= 9----------------------------------
par(mfrow=c(2,2),mar = c(2,3,2,2))
# Plot the fitted values
plot(fitted(fit4, "full"), main = "Fitted Values")

# Check Residuals
plot(residuals(fit4, "raw"), main = "Raw Residuals")

plot(residuals(fit4, "pearson"), main = "Pearson Residuals")

plot(residuals(fit4, "deviance"), main = "Deviance Residuals")

## ---- fig.height= 4, fig.width= 4----------------------------------------
head(AnxStrData, 8)

plot(density(AnxStrData$Anxiety), main = "Anxiety and Stress")
lines(density(AnxStrData$Stress), lty = 2)

## ----fit2----------------------------------------------------------------
# Fit the null model
fit_null <- cdfquantreg(Anxiety ~ 1 | 1, fd, sd, data = AnxStrData)

# Fit the location model
fit1 <- cdfquantreg(Anxiety ~ Stress | 1, fd, sd, data = AnxStrData)

# Fit the full model
fit2 <- cdfquantreg(Anxiety ~ Stress | Stress, fd, sd, data = AnxStrData)

anova(fit_null,fit1, fit2)

summary(fit2)

## ----plotfit2, fig.height= 4, fig.width= 4-------------------------------
# Compare the empirical distribution and the fitted values distribution
plot(fit2)

# Plot the fitted values
plot(fitted(fit2, "full"))

# Check Residuals
plot(residuals(fit2, "raw"))


