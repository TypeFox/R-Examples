###hop:3-9
summary(eband.model)
4.554 + c(-1,1) * 1.543 * qt(0.975,df=5)          # CI by hand
confint(eband.model,"stretch")                   # CI using confint()
