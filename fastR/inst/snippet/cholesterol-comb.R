model1 <- lm(response~trt,cholesterol)
cholesterol$x1 <- cholesterol$trt == "drugD"
cholesterol$x2 <- cholesterol$trt == "drugE"
model2 <- lm(response~ 1 + x1 + x2 , cholesterol)
anova(model1,model2)
