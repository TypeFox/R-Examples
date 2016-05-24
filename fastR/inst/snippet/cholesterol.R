data(cholesterol,package="multcomp")
chol.model <- lm(response~trt,cholesterol)
plot(chol.model)              # diagnostic plots
###hop:3-9
summary(chol.model)
anova(chol.model)
