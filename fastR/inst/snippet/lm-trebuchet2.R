treb.model <- lm(distance~projectileWt,data=trebuchet2)
coef(treb.model)
plot1 <- xyplot(distance~projectileWt,data=trebuchet2, type=c('p','r'))
