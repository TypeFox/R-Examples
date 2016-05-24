require(alr3); data(rat)
rat.lm <- lm(y~BodyWt*LiverWt,rat)
summary(rat.lm)
