airp.lm1 <- lm(pollution~location,airpollution)
# specify contrasts by giving the coefficients
contr <- rbind(
	c(0,1,0),
	c(0,0.5,-1))
# we can give our contrasts custom names if we like
contr1 <- rbind(
	"hill - plains" = c(0,1,0),
	"suburb - urban" = c(0,0.5,-1))
###hop:3-6
summary(glht(airp.lm1,contr1))
