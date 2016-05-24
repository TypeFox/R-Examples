poison.lm <- lm(Time~factor(Poison) * factor(Treatment),poison)
anova(poison.lm)
