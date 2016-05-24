require(faraway); data(rabbit,package="faraway")
rabbit.lm1 <- lm(gain~block+treat,rabbit)
rabbit.lm2 <- lm(gain~treat+block,rabbit)
rabbit.lm3 <- lm(gain~block,rabbit)
rabbit.lm4 <- lm(gain~treat,rabbit)
anova(rabbit.lm1)
anova(rabbit.lm2)
anova(rabbit.lm1,rabbit.lm3)
anova(rabbit.lm2,rabbit.lm4)
# width of Tukey HSD intervals:
2 * qtukey(0.95,6,15) / sqrt(5) * summary(rabbit.lm1)$sigma
# TukeyHSD() can automate this:
TukeyHSD(aov(gain~block+treat,rabbit), "treat")
