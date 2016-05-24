summary( abs(mle.count-lambda) < abs(mle.measure - lambda) ~ size+lambda,
    data=results,
    method="cross",
    fun = function(x) { round(mean(x),3) } 
    )

plot1a <- bwplot( (mle - lambda)/lambda ~ method | factor(lambda),
    subset=size==50,
    data=results2, 
    groups=method,
    auto.key=T)

plot2 <- xyplot( (mle - lambda)/lambda ~ area | factor(lambda)*factor(size), 
    data=results2, 
    groups=method, 
    scales=list(relation="free"),
    panel = function(x,y,...) {
        panel.abline(h=0,lwd=2,col="red")
        panel.xyplot(x,y,...)
        },
    auto.key=T,
    layout=c(3,2),
    )

plot3 <- xyplot( mle - lambda ~ plants | factor(lambda)*factor(size), 
    data=results2, 
    groups=method, 
    scales=list(relation="free"),
    panel = function(x,y,...) {
        panel.abline(h=0,lwd=2,col="red")
        panel.xyplot(x,y,...)
        },
    auto.key=T,
    layout=c(3,2),
    )
