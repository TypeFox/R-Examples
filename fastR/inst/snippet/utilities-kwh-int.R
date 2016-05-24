ut.lm3 <- lm(thermsPerDay ~ temp * kwhpday, ut)
###hop:3-9
summary(ut.lm3)
ut.plot3 <- xplot(ut.lm3)
