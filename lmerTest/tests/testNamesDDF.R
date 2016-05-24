require(lmerTest)


m <- lmer(Informed.liking ~ Gender+Information+Product +(1|Consumer), data=ham)


anova(m)
anova(m, ddf="sat")
anova(m, ddf="SA")
anova(m, ddf="kenw")
anova(m, ddf="lme4")
