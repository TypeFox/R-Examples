choose(90,4)/choose(100,4)  # prob only good ones selected
1 - choose(90,4)/choose(100,4)  # lot is rejected
f <- function(x) { 1 - choose(100-x,4)/choose(100,4) }
myplot <- xyplot(sapply(10:100,f) ~ 10:100, type='l',
        xlab="number of defective parts",
        ylab="probability of rejecting",
        lwd=2,
        col="navy")
