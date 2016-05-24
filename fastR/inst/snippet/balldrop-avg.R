balldropavg <- aggregate(balldrop$time, by=list(balldrop$height), mean)
names(balldropavg) <- c('height','time')
ball.modelA <- lm(time ~ sqrt(height),balldropavg)
###hop:3-9
summary(ball.modelA)
ball.plotA <- xyplot(time~height,balldropavg,
                panel=panel.lm,model=ball.modelA)
ball.residplotA <- xplot(ball.modelA,w=1)
