ball.model <- lm(time~height,balldrop)
###hop:3-9
summary(ball.model)
ball.plot <- xyplot(time~height,balldrop,type=c('p','r'))
ball.residplot <- xplot(ball.model,w=1)
