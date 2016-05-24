model1 <- lm(score~noise+group,mathnoise)
anova(model1)
summary(score~group,data=mathnoise,fun=favstats)
