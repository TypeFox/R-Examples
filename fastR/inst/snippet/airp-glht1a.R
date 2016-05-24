# these look nicer if we parameterize differently in the model
airp.lm2 <- lm(pollution~-1 + location,airpollution)
contr2 <- rbind(
	"hill - plains" = c(1,-1,0),
	"suburb - urban" = c(1,1,-2))
###hop:3-6
summary(glht(airp.lm2,contr2))
