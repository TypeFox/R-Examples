# remove first few observations because of bad meter read
ut2 <- subset( utilities, subset=(year > 2000 | month > 6) )
ut2.lm <- lm( thermsPerDay ~ month + I(month^2), ut2 )
###hop:3-9
summary(ut2.lm)
ut2.plot1 <- xyplot(thermsPerDay ~ month, data=ut2,
                      panel=panel.lm, model=ut2.lm)
ut2.plot2 <- xplot(ut2.lm, w=1:2)
