ball.modelT <- lm(time ~ sqrt(height),balldrop)
###hop:3-9
summary(ball.modelT)
ball.plotT <- xyplot(time~height,balldrop,
                panel=panel.lm,model=ball.modelT)
ball.residplotT <- xplot(ball.modelT,w=1)
