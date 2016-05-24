rat.lm1 <- lm(consumption~flavor,ratpoison)
anova(rat.lm1)
summary(rat.lm)$sigma
summary(rat.lm1)$sigma
(summary(rat.lm1)$sigma/summary(rat.lm)$sigma)^2
###hop:3-9
summary(rat.lm)
###hop:3-9
summary(rat.lm1)
