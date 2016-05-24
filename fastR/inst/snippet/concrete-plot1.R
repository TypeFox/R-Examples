plot1 <- xplot(concrete.lm0,which=2)
plot2 <- xplot(concrete.lm0,which=3)
plot3 <- xplot(concrete.lm0,which=5,add.smooth=FALSE)
plot4 <- xyplot(resid(concrete.lm1) ~ fitted(concrete.lm1),
               main = "residuals vs fits",
               ylab = "residuals",
               xlab = "fitted values",
               sub  = "lm(strength ~ limestone + water)" )
