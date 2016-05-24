summary(lm(weight~groove,tirewear))
plot <- xyplot(weight~groove,tirewear, type=c('p','r'))
