# using mcp() to help build the contrasts
airp.lm3 <- lm(pollution~location,airpollution)
contr3 <- mcp(location = rbind(
	"hill - plains" = c(1,-1,0),
	"suburb - urban" = c(1,1,-2)
	))
###hop:3-9
summary(glht(airp.lm3,contr3))
