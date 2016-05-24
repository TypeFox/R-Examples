chol.glht <- confint(glht(chol.model,mcp(trt="Tukey")))
###hop:3-9
summary(chol.glht)
plot(confint(chol.glht))
