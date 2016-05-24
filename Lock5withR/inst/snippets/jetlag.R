jetlag.model <- lm( shift ~ treatment, JetLagKnees )
anova(jetlag.model)
summary(jetlag.model)
plot(jetlag.model, w = 1:2)

