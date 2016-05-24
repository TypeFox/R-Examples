# subset the data:
#   * remove first few months where there appears to have been a bad
#         meter reading (year == 2000 & month <= 6)
#   * remove months where there is little need to heat (temp > 60)
ut <- subset(utilities2, subset=(year > 2000 | month > 6) & temp <= 60)
ut.lm1 <- lm(thermsPerDay ~ temp, ut)
###hop:3-9
summary(ut.lm1)
ut.plot <- xplot(ut.lm1)
