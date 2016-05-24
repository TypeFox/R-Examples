pow <- function(effect) {
    power.t.test(delta=effect,n=50)$power
}
effect = seq(0,2,by=0.05)
plot <- xyplot(pow(effect) ~ effect, type='l',
    ylab="power", xlab="effect size",
    main="Power of a 2-sample test (n=50)")
