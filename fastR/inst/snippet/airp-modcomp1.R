# convert location to a numeric variable for convenience
airp <- airpollution; airp$loc <- as.numeric(airp$location); airp
model <- lm(pollution~location, airp)
model2 <- lm(pollution~  1 + (loc==3), airp)
anova(model2,model)
