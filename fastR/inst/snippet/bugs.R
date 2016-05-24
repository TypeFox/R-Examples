model <- aov(sqrt(NumTrap)~Color,bugs)
TukeyHSD(model)
model <- lm(sqrt(NumTrap)~Color,bugs)
summary(glht(model,mcp(Color="Tukey")))
