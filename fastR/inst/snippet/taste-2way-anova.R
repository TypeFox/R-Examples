taste.lm <- lm(score~scr*liq,data=tastetest)
anova(taste.lm)
