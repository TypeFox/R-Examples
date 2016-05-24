poison.lm2 <- lm(1/Time~factor(Poison) * factor(Treatment),poison)
xplot(poison.lm2,w=c(4,2))
