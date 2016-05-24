require(multcomp)
coag.glht <- glht(coag.model,mcp(diet="Tukey"))
summary(coag.glht)  
plot(TukeyHSD(coag.aov));  plot(confint(coag.glht))   # plots
