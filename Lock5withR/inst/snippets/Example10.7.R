Mod0 <- lm(Price ~ 1, data = InkjetPrinters)
Mod1 <- lm(Price ~ PhotoTime + CostColor, data = InkjetPrinters)
msummary(Mod1)
anova(Mod0, Mod1)

