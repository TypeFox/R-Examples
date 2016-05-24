require(faraway); data(eggprod,package="faraway")
eggprod.lm <- lm(eggs~block+treat,eggprod)
anova(eggprod.lm)
###hop:3-9
summary(eggprod.lm)
coef(summary(eggprod.lm))
eggprod.lm1way <- lm(eggs~treat,eggprod)
anova(eggprod.lm1way)
# width of Tukey HSD intervals:
2 * qtukey(0.95,3,6) / sqrt(4) * summary(eggprod.lm)$sigma
# TukeyHSD() can automate this:
TukeyHSD(aov(eggs~block+treat,eggprod), "treat")
