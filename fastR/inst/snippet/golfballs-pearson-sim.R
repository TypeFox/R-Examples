E <- rep(486/4,4)
rgolfballs <- rmultinom(n=10000,size=486,prob=rep(0.25,4))
chisqstat <- function(x) { sum((x-E)^2 / E) }
plot3 <- statTally(golfballs, rgolfballs, chisqstat, 
                   xlab=expression(X^2))
