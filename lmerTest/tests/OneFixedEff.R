require(lmerTest)

m <- lmer(Coloursaturation ~ TVset+
            (1|Assessor), data=TVbo)

step(m)


m2 <- lmer(Coloursaturation ~ 
            (1|Assessor), data=TVbo)

anova(m2)

summary(m2)
