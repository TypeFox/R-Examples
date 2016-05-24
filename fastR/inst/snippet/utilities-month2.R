ut2$monthShifted <- ( ut2$month -2 ) %% 12
ut2.lm2 <- lm(thermsPerDay ~ monthShifted + I(monthShifted^2), ut2)
###hop:3-9
summary(ut2.lm2)
ut2.plot3 <- xyplot(thermsPerDay ~ monthShifted, data=ut2,
                      panel=panel.lm, model=ut2.lm2)
ut2.plot4 <- xplot(ut2.lm2, w=1:2)
