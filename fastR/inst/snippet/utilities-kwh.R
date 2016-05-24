ut.lm2 <- lm(thermsPerDay ~ temp + kwhpday, ut)
###hop:3-9
summary(ut.lm2)
ut.plot2 <- xplot(ut.lm2)
