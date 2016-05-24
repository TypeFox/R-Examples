summary(score~type,data=tastetest,
	fun=function(x){cbind(mean=mean(x),sd=sd(x))})
taste.xy <- xyplot(score~type,data=tastetest)
taste.lm <- lm(score~type,data=tastetest)
anova(taste.lm)
###hop:3-12
taste.cint <- confint(glht(taste.lm,mcp(type="Tukey"))); taste.cint
